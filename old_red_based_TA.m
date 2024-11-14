clc
clear
close all

%% YALMIP

addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master');
addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master/extras');
addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master/solvers');
addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master/modules');
addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master/modules/parametric');
addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master/modules/moment');
addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master/modules/global');
addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master/modules/sos');
addpath('/home/lorenzo/Documents/MATLAB/YALMIP-master/operators');

addpath('/home/lorenzo/Documents/MATLAB/sedumi-master');

%% Adjacency matrix

ad = [0, 0, 0, 0, 1, 0, 0;
      0, 0, 1, 0, 1, 0, 0;
      0, 1, 0, 1, 0, 0, 1;
      0, 0, 1, 0, 0, 1, 0;
      1, 1, 0, 0, 0, 0, 0;
      0, 0, 0, 1, 0, 0, 0;
      0, 0, 1, 0, 0, 0, 0];

%% Optimization parameters

l = 5;
C = 50;

%% Constraints on how many robot are needed for a particular task

% Task 1 --> coverage

n1_min = 2; 
n1_max = 3;

% Task 2 --> point in the plane

n2_min = 2;
n2_max = 3;

% Task 3 --> hovering position

n3_min = 1;
n3_max = 1;

%% Mappings 

nt = 3; % # of tasks
nc = 3; % # of capabilities
nf = 4; % # of features
nr = 7; % # of robots

% Capability-to-Task Mapping (nt x nc)

T = [1 0 0;
     0 1 0;
     0 0 1];

% Feature-to-Capability Mapping (we have three capabilities)

H1 = [0.5 0.5 0 0;
     0 0 0 1];

H2 = [0.5 0 0.5 0];

H3 = [0 0 0 1];

% Robot-to-Feature Mapping

A_old = [1 1 1 1 1 0 0;
     1 1 0 1 0 0 0;
     0 0 1 0 1 0 0;
     0 0 0 0 0 1 1];

% 1-hop sensor grouping

A = sensor_grouping(A_old,ad,1);

% Robot-to-Capability Mapping

F = Rob2Cap_Matrix(H1,H2,H3,A);

% Specialization matrix (nt x nt)

S1 = specialization(T,F(:,1),nt);
S2 = specialization(T,F(:,2),nt);
S3 = specialization(T,F(:,3),nt);
S4 = specialization(T,F(:,4),nt);
S5 = specialization(T,F(:,5),nt);
S6 = specialization(T,F(:,6),nt);
S7 = specialization(T,F(:,7),nt);

% Redundancy matrices (nr*nr)

% Matrix 1
R1 = [0   0   0   0   2/3 0   0;
      0   0   2.0 0   2/5 0   0;
      0   2.0   0   2.0 0 0   2.0;
      0   0   2.0 0   0   2.0   0;
      2/3 2/5 0   0   0   0   0;
      0   0   0   2.0   0   0   0;
      0   0   2.0 0   0   0   0];

% R1 = R1 + 5*eye(size(R1));

% Matrix 2
R2 = [0   0   0   0   2.0 0   0;
      0   0   2.0   0   2.0 0   0;
      0   2.0   0   2/3 0   0   2/5;
      0   0   2/3 0   0   2.0   0;
      2.0 2.0 0   0   0   0   0;
      0   0   0   2.0   0   0   0;
      0   0   2/5 0   0   0   0];

% R2 = R2 + 5*eye(size(R2));

% Matrix 3
R3 = [zeros(1,7);
      zeros(1,7);
      0 0   0   0 0 0 2.0;
      0 0 0 0 0 2.0 0;
      zeros(1,7);
      0 0   0   2.0   0 0 0;
      0 0   2.0 0   0 0 0];
% 
% R3 = R3 + 5*eye(size(R3));
% R3(6,6) = 2/5;

% % Matrix 1
% R1 = [0   0.7 0   0   0.7 0 0;
%       0.7 0   0   0   0   0 0;
%       0   0   0   0.3 0   0 0.2;
%       0   0   0.3 0   0   0 0;
%       0.7 0   0   0   0   0 0;
%       0   0   0   0   0   0 0;
%       0   0   0.2 0   0   0 0];
% 
% % Matrix 2
% R2 = [0   0.5 0   0   0.2 0 0;
%       0.5 0   0   0   0   0 0;
%       zeros(1,7);
%       zeros(1,7);
%       0.2 0   0   0   0   0 0;
%       zeros(1,7);
%       zeros(1,7)];
% 
% % Matrix 3
% R3 = [zeros(1,7);
%       zeros(1,7);
%       0 0   0   0.9 0 0 0.3;
%       0 0   0.9 0   0 0 0;
%       zeros(1,7);
%       0 0   0   0   0 0 0;
%       0 0   0.3 0   0 0 0];

sum(R1,2)

sum(R2,2)

sum(R3,2)

red = blkdiag(R1,R2,R3);

%% Constraints

% Matrices for alpha constraint (2b)
    
C3 = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1];

c3 = ones(nr,1);

% Matrices for constraint (2c)

F_t1 = constraint_11_e(nc,nr,nt,F,1);
F_t2 = constraint_11_e(nc,nr,nt,F,2);
F_t3 = constraint_11_e(nc,nr,nt,F,3);

C4 = -[F_t1;F_t2;F_t3];

c4 = [-T(1,:)';-T(2,:)';-T(3,:)'];

% Matrices for number of robots constraint (2d)

C5 = [1 0 0, 1 0 0, 1 0 0, 1 0 0, 1 0 0, 1 0 0, 1 0 0;
     -1 0 0, -1 0 0, -1 0 0, -1 0 0, -1 0 0, -1 0 0, -1 0 0;
      0 1 0, 0 1 0, 0 1 0, 0 1 0, 0 1 0, 0 1 0, 0 1 0;
      0 -1 0, 0 -1 0, 0 -1 0, 0 -1 0, 0 -1 0, 0 -1 0, 0 -1 0;
      0 0 1, 0 0 1, 0 0 1, 0 0 1, 0 0 1, 0 0 1, 0 0 1;
      0 0 -1, 0 0 -1, 0 0 -1, 0 0 -1, 0 0 -1, 0 0 -1, 0 0 -1
     ];

c5 = [n1_max; -n1_min; n2_max; -n2_min; n3_max; -n3_min];

%% Cost function

% Pi matrix

P1 = eye(nt)-S1*pinv(S1);
P2 = eye(nt)-S2*pinv(S2);
P3 = eye(nt)-S3*pinv(S3);
P4 = eye(nt)-S4*pinv(S4);
P5 = eye(nt)-S5*pinv(S5);
P6 = eye(nt)-S6*pinv(S6);
P7 = eye(nt)-S7*pinv(S7);

P1_bar = P1'*P1;
P2_bar = P2'*P2;
P3_bar = P3'*P3;
P4_bar = P4'*P4;
P5_bar = P5'*P5;
P6_bar = P6'*P6;
P7_bar = P7'*P7;

p = blkdiag(P1_bar,P2_bar,P3_bar,P4_bar,P5_bar,P6_bar,P7_bar);

% Bad allocation term

H1 = C*p;

% alpha matrix

alpha = binvar(nt*nr,1);

%% Task Allocation

Z1_task = Task_All(C3,C4,C5,c3,c4,c5,alpha,H1,R1,R2,R3);

Alpha1(1,1:nt) = Z1_task(1:3); Alpha1(2,1:nt) = Z1_task(4:6);
Alpha1(3,1:nt) = Z1_task(7:9); Alpha1(4,1:nt) = Z1_task(10:12);
Alpha1(5,1:nt) = Z1_task(13:15); Alpha1(6,1:nt) = Z1_task(16:18);
Alpha1(7,1:nt) = Z1_task(19:21);

Alpha1

%% Functions

function z = Task_All(C3,C4,C5,c3,c4,c5,alpha,H1,R1,R2,R3)

ops = sdpsettings('verbose',0);

% Inequality constraints matrices

A3 = C3*alpha <= c3;
A4 = C4*alpha <= c4;
A5 = C5*alpha <= c5;

% YALMIP

F = [A3; A4; A5];

obj1 = alpha'*H1*alpha;

first_row = [alpha(1); alpha(4); alpha(7); alpha(10); alpha(13); alpha(16); alpha(19)];
second_row = [alpha(2); alpha(5); alpha(8); alpha(11); alpha(14); alpha(17); alpha(20)];
third_row = [alpha(3); alpha(6); alpha(9); alpha(12); alpha(15); alpha(18); alpha(21)];

obj2 = first_row'*(R1')*R1*first_row + second_row'*(R2')*R2*second_row + third_row'*(R3')*R3*third_row;

obj = 50*obj1 + 25*obj2;

optimize(F,obj,ops);

z = value(alpha);


end

function y = max_kron(n,x)
%KRON Summary of this function goes here
%   Detailed explanation goes here

[row,col] = size(x);

temp = zeros(row,col);

y = zeros(1,col);

for i = 1:row    
    for j = 1:col       
        if ( x(i,j)==n )            
            temp(i,j) = 1;            
        end        
    end
end

for j = 1:col  
    max = 0;
    for i =1:row
        if ( temp(i,j) >= max )
            max = temp(i,j);
        end
    end
    y(j) = max;
end

end


function F = Rob2Cap_Matrix(H1,H2,H3,A)

f1 = H1*A;
f2 = H2*A;
f3 = H3*A;

F1 = max_kron(1,f1);
F2 = max_kron(1,f2);
F3 = max_kron(1,f3);

F = [F1; F2; F3];

end

function F_t = constraint_11_e(nc,nr,nt,F,index)

F_t = zeros(nc, nt*nr);

j = 1;

for i = 1:nc
    
    h = index;
    
    for k = 1:nr
        
        F_t(i,h) = F(j,k);
        h = h+nt;

    end

    j = j+1;

end

end

function S = specialization(T,F,nt)
%SPECIALIZATION Summary of this function goes here
%   Detailed explanation goes here

S = zeros(nt,nt);

for i=1:nt
    
    if ( T(i,:)*F > 0 )
        S(i,i) = 1;
    end

end
end

