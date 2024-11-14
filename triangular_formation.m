% Execution of the triangular formation

delta_x_t1 = [0, -2, 2]';
delta_y_t1 = [0, 2, 2]';

ad_t1 = ad(cc{idx1},cc{idx1});

D_t1 = diag(sum(ad_t1,1));

Lt1 = D_t1 - ad_t1;

u_max = [u_max_x(cc{idx1});u_max_y(cc{idx1})];

[x_nom_t1,y_nom_t1,yc_t1,tspan_t1] = consensus(Lt1,x0_cc1,y0_cc1,dt,Tf,delta_x_t1,delta_y_t1);

% figure()
% hold on; grid on;
% plot(x_nom_t1(:,1),y_nom_t1(:,1),'o');
% plot(x_nom_t1(:,end),y_nom_t1(:,end),'*');
% plot(x_nom_t1(1,:),y_nom_t1(1,:));
% plot(x_nom_t1(2,:),y_nom_t1(2,:));
% plot(x_nom_t1(3,:),y_nom_t1(3,:));

[x_sat_t1,y_sat_t1,yc_sat_t1,tspan_sat_t1] = consensus_sat(Lt1,x0_cc1,y0_cc1,dt,Tf,delta_x_t1,delta_y_t1,u_max);

% figure()
% hold on; grid on;
% plot(x_sat_t1(:,1),y_sat_t1(:,1),'o');
% plot(x_sat_t1(:,end),y_sat_t1(:,end),'*');
% plot(x_sat_t1(1,:),y_sat_t1(1,:));
% plot(x_sat_t1(2,:),y_sat_t1(2,:));
% plot(x_sat_t1(3,:),y_sat_t1(3,:));

% AGENT DYNAMICS

% number of robots
n_t1 = max(size(Lt1));

% dynamical matrices
A = zeros(2*n_t1,2*n_t1);
B = eye(2*n_t1);
C = blkdiag(Lt1,Lt1);

% ALLOCATOR
s = tf('s');

% annihilator
N = C;
N_perp = null(C);
g = s+2;
W_an = zpk(minreal((N_perp*(s+1))/g));

ss_An = ss(W_an);

A_an = ss_An.A;
B_an = ss_An.B;
C_an = ss_An.C;
D_an = ss_An.D;

W_an0 = dcgain(W_an);
x_an0 = zeros(max(size(A_an)),1);

% optimizer 

alpha = 50;
C_opt = eye(2);
x_opt0 = zeros(2,1);

% threshold

th = 0.05;

tf = 100;
f = 0.2;
delta = [delta_x_t1;delta_y_t1];
x0 = x0_cc1;
y0 = y0_cc1;
num = max(size(x0));
sig_vec = ones(num,1);

out_dead_t1 = sim('consensus_dead');

t_dead_t1 = out_dead_t1.t;
dim_t1 = max(size(t_dead_t1));
x_dead_t1 = zeros(num,dim_t1);
y_dead_t1 = zeros(num,dim_t1);
x_dead_t1(:,:) = out_dead_t1.y(1:num,:,:);
y_dead_t1(:,:) = out_dead_t1.y(num+1:end,:,:);
u_sat_dead_t1 = zeros(2*num,dim_t1);
u_sat_dead_t1(:,:) = out_dead_t1.u(:,:,:);

% figure()
% hold on; grid on;
% plot(x_dead_t1(:,1),y_dead_t1(:,1),'o');
% plot(x_dead_t1(:,end),y_dead_t1(:,end),'*');
% plot(x_dead_t1(1,:),y_dead_t1(1,:));
% plot(x_dead_t1(2,:),y_dead_t1(2,:));
% plot(x_dead_t1(3,:),y_dead_t1(3,:));
