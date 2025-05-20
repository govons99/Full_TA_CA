clc
clear
close all

addpath('/home/lorenzo/Documents/PhD/Second_year/Full_TA_CA/simulation')

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

ad = ad+eye(10);

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

%% New adjacency matrix

ad_new = zeros(np,np);

ad_new(1:10,1:10) = ad;
ad_new(1:10,11:end) = M(:,11:end);
ad_new(11:end,1:10) = (M(:,11:end))';

for i=11:np
    for j=11:np
        if i~=j
            if intersect(cc{i},cc{j})
                ad_new(i,j) = 1;
            end
        end
    end
end

ad_new = ad_new + eye(np);

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

%% Redundancy matrices (nr*nr)

scenario = 1;

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

%% Bid vectors

w1 = 1;
w2 = 2;

N = 10;

r1 = diag(zeros(N,1));
r2 = diag(zeros(N,1));
r3 = diag(zeros(N,1));

r1(5,5) = 1;
r1(8,8) = 1;
r1(9,9) = 1;

r2(2,2) = 2;
r2(3,3) = 2;
r2(4,4) = 2;
r2(6,6) = 2;

r3(1,1) = 1;
r3(7,7) = 1;
r3(10,10) = 1;

r = {r1,r2,r3};

c = zeros(N,nt);

for i=1:N
    for m=1:nt
        if 1
            c(i,m) = w1 * 1 + w2 * (r{m}(i,i));
        end
    end
end

%% CAAB

x_ini = zeros(N,nt);
y_ini = zeros(N,nt);
J_ini = zeros(N,1);

x_ph1 = zeros(N,nt);
y_ph1 = zeros(N,nt);

x_ph2 = zeros(N,nt);
y_ph2 = zeros(N,nt);

J = zeros(N,1);

tic;
for t=1:5
    disp('#########################')
    % phase 1
    for i=1:N

        [xi_1,yi_1,Ji] = phase1(c(i,:)',x_ini(i,:),y_ini(i,:),J_ini(i));

        x_ph1(i,:) = xi_1;
        y_ph1(i,:) = yi_1;
        J(i) = Ji;

        J_ini(i) = Ji;

    end
    
    % phase 2
    for i=1:N

        mask = logical(ad(i,:));         % your selection mask
        result = zeros(size(y_ph1));         % initialize with zeros
        result(mask, :) = y_ph1(mask, :);

        [xi_2,yi_2] = phase2(x_ph1(i,:),y_ph1(i,:),y_ph1(logical(ad(i,:)),:),J(i),i);

        x_ph2(i,:) = xi_2;
        y_ph2(i,:) = yi_2;

        x_ini(i,:) = xi_2;
        y_ini(i,:) = yi_2;

    end
    x_ph2

end
toc
