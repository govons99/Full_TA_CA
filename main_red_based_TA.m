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

addpath('/opt/gurobi1103/linux64/matlab')
savepath

%% Communication graph and clusterization

% adjacency matrix

ad = [0 1 0 0 0 0 0 0 0 1;
      1 0 1 1 0 0 0 0 0 0;
      0 1 0 0 0 1 0 0 1 0;
      0 1 0 0 0 1 1 0 1 0;
      0 0 0 0 0 0 0 1 0 1;
      0 0 1 1 0 0 0 0 0 0;
      0 0 0 1 0 0 0 1 0 1;
      0 0 0 0 1 0 1 0 1 0;
      0 0 1 1 0 0 0 1 0 0;
      1 0 0 0 1 0 1 0 0 0];

D = diag(sum(ad,1));

L = D - ad;

% number of agents to put in a cluster

N = 4;

% connected components

cc = {};
p = [];

for i=1:N
    comp = connected_components(ad,i);
    dim = max(size(cc));
    count = 1;
    for j = dim+1:dim+max(size(comp))
        cc{j} = comp{count};
        p(j) = max(size(cc{j}));
        count = count+1;
    end
end

ix = 1;
b = zeros(dim*4,1);
s = 's';
for i=5:10
    for j=1:max(size(cc))
        if ismember(i,cc{j})
            h(ix) = i;
            s = num2str(cc{j});
            ix = ix+1;
        end
    end
end

% parameters

na = max(size(ad));
np = max(size(cc));

% Robot-to-Partition Mapping

M = zeros(na,np);
M(1:na,1:na) = eye(na);

for i = 1:na
    for j = na+1:np
          if any(cc{j}==i)
              M(i,j) = 1;
          end
    end
end


%% Optimization parameters

l = 5;
c = 50;

%% Mappings 

nt = 3; % # of tasks
nc = 3; % # of capabilities
nf = 4; % # of features
nr = np; % # of robots

% Capability-to-Task Mapping (nt x nc)

T = eye(nt);

% Feature-to-Capability Mapping (we have three capabilities)

H1 = [0.5 0.5 0 0;
      0 0 0 1];

H2 = [0.5 0 0.5 0];

H3 = [0 0 0 1];

H = {H1,H2,H3};

% Robot-to-Feature Mapping

A_old = [0 1 0 1 1 0 1 1 1 1
         0 0 0 0 1 0 1 0 1 0;
         0 1 0 1 0 0 0 1 0 1;
         1 0 1 0 0 1 0 0 0 0];

% Sensor grouping

A = [A_old A_old*M(:,na+1:end)];

for i=1:nf
    for j=1:np
        if A(i,j)
            A(i,j) = 1;
        end
    end
end

% Robot-to-Capability Mapping

F = Rob2Cap_Matrix(H,A,nc,nr);

% Specialization matrix (nt x nt)

S = {};

for i=1:nr
    S{i} = specialization(T,F(:,i),nt);
end

% Pi matrix

Pi = {};

for i=1:nr
    Pi{i} = eye(nt)-S{i}*pinv(S{i});
end

%% Integer variables

alpha = binvar(nt,nr);
delta = binvar(nr,1);

%% Task Allocation

%scenario = input('Choose a scenario (is an integer number between 1 and 3): ');

scenario = 1;

% Redundancy matrices (nr*nr)

R1 = diag(ones(nr,1));
R2 = diag(ones(nr,1));
R3 = diag(ones(nr,1));

switch scenario
    case 1
        R1(48,48) = 0.2;
        R2(66,66) = 0.4;
        R3(29,29) = 0.6;
    case 2
        R1(28,28) = 0.2;
        R2(73,73) = 0.4;
        R3(38,38) = 0.6;
    case 3
        R1(46,46) = 0.2;
        R2(55,55) = 0.4;
        R3(41,41) = 0.6;
    otherwise
        disp('wrong number inserted')
        disp('no redundancy case')
end

R = {R1,R2,R3};

% Constraints on how many robot are needed for a particular task

% Task 1 --> triangle formation

n1_min = 3; 
n1_max = 3;

% Task 2 --> square formation

n2_min = 4;
n2_max = 4;

% Task 3 --> rendezvous

n3_min = 3;
n3_max = 3;

% min max robot vectors

n_min = [n1_min,n2_min,n3_min];
n_max = [n1_max,n2_max,n3_max];

% Allocation

[Z1_task, Delta, elapsed] = Task_All_new(nr,Pi,alpha,nt,R,delta,F,T,n_min,n_max,p,c,l,M);

disp('Task t1 has been assigned to:')
idx1 = find(Z1_task(1,:));
cc{idx1}

disp('Task t2 has been assigned to:')
idx2 = find(Z1_task(2,:));
cc{idx2}

disp('Task t3 has been assigned to:')
idx3 = find(Z1_task(3,:));
cc{idx3}

disp(['The allocation took: ',num2str(elapsed),' seconds'])

%% Task execution

% time interval

dt = 0.01;

Tf = 100;

% initial condition

x0 = zeros(na,1);
y0 = zeros(na,1);

for i=1:na
    x0(i) = 10*cos((2*pi)/na*i);
    y0(i) = 10*sin((2*pi)/na*i);
end

% saturation

u_max_x = [1,0.5,1,0.4,1.5,0.5,1,0.5,1,0.4]';
u_max_y = [1,0.5,1,0.4,1.5,0.5,1,0.5,1,0.4]';

% initial condition clusters

x0_cc1 = x0(cc{idx1});
y0_cc1 = y0(cc{idx1});

x0_cc2 = x0(cc{idx2});
y0_cc2 = y0(cc{idx2});

x0_cc3 = x0(cc{idx3});
y0_cc3 = y0(cc{idx3});

%% execution of the triangular formation

triangular_formation

%% execution of the squared formation

squared_formation_dec

%% execution of the rendez-vous

rendez_vous

%% plot

plots

%% txt files

%txt_writer
