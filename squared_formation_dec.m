% execution of the squared formation

delta_x_t2 = [0, 0, 2, 2]';
delta_y_t2 = [0, 2, 0, 2]';

ad_t2 = ad(cc{idx2},cc{idx2});

D_t2 = diag(sum(ad_t2,1));

Lt2 = D_t2 - ad_t2;

u_max = [u_max_x(cc{idx2});u_max_y(cc{idx2})];

%% nominal case

[x_nom_t2,y_nom_t2,yc_t2,tspan_t2] = consensus(Lt2,x0_cc2,y0_cc2,dt,Tf,delta_x_t2,delta_y_t2);

% figure()
% hold on; grid on;
% plot(x_nom_t2(:,1),y_nom_t2(:,1),'o');
% plot(x_nom_t2(:,end),y_nom_t2(:,end),'*');
% plot(x_nom_t2(1,:),y_nom_t2(1,:));
% plot(x_nom_t2(2,:),y_nom_t2(2,:));
% plot(x_nom_t2(3,:),y_nom_t2(3,:));
% plot(x_nom_t2(4,:),y_nom_t2(4,:));

%% saturation case

[x_sat_t2,y_sat_t2,yc_sat_t2,tspan_sat_t2] = consensus_sat(Lt2,x0_cc2,y0_cc2,dt,Tf,delta_x_t2,delta_y_t2,u_max);

% figure()
% hold on; grid on;
% plot(x_sat_t2(:,1),y_sat_t2(:,1),'o');
% plot(x_sat_t2(:,end),y_sat_t2(:,end),'*');
% plot(x_sat_t2(1,:),y_sat_t2(1,:));
% plot(x_sat_t2(2,:),y_sat_t2(2,:));
% plot(x_sat_t2(3,:),y_sat_t2(3,:));
% plot(x_sat_t2(4,:),y_sat_t2(4,:));

%% decentralized steady-state optimization

% neighbours

idx_n2 = [3 2 4];
idx_n3 = [2 3 6];
idx_n6 = [3 6 4];
idx_n4 = [6 4 2];

% index in the cluster

Idx2 = [2 1 3];
Idx3 = [1 2 4];
Idx6 = [2 4 3];
Idx4 = [4 3 1];

% laplacians

ad_n2 = ad(idx_n2,idx_n2);
D_n2 = diag(sum(ad_n2,1));
L_n2 = D_n2 - ad_n2;

ad_n3 = ad(idx_n3,idx_n3);
D_n3 = diag(sum(ad_n3,1));
L_n3 = D_n3 - ad_n3;

ad_n6 = ad(idx_n6,idx_n6);
D_n6 = diag(sum(ad_n6,1));
L_n6 = D_n6 - ad_n6;

ad_n4 = ad(idx_n4,idx_n4);
D_n4 = diag(sum(ad_n4,1));
L_n4 = D_n4 - ad_n4;

% number of robots
n_t2 = max(size(L_n2));

% dynamical matrices
A = zeros(2*n_t2,2*n_t2);
B = eye(2*n_t2);

C_n2 = blkdiag(L_n2,L_n2);
C_n3 = blkdiag(L_n3,L_n3);
C_n6 = blkdiag(L_n6,L_n6);
C_n4 = blkdiag(L_n4,L_n4);

% annihilator neighbour 2
N = C_n2;
N_perp = null(N);
g = s+2;
W_an = zpk(minreal((N_perp*(s+1))/g));

ss_An = ss(W_an);

A_an_n2 = ss_An.A;
B_an_n2 = ss_An.B;
C_an_n2 = ss_An.C;
D_an_n2 = ss_An.D;

W_an0_n2 = dcgain(W_an);
x_an0 = zeros(max(size(A_an_n2)),1);

% annihilator neighbour 3
N = C_n3;
N_perp = null(N);
g = s+2;
W_an = zpk(minreal((N_perp*(s+1))/g));

ss_An = ss(W_an);

A_an_n3 = ss_An.A;
B_an_n3 = ss_An.B;
C_an_n3 = ss_An.C;
D_an_n3 = ss_An.D;

W_an0_n3 = dcgain(W_an);

% annihilator neighbour 6
N = C_n6;
N_perp = null(N);
g = s+2;
W_an = zpk(minreal((N_perp*(s+1))/g));

ss_An = ss(W_an);

A_an_n6 = ss_An.A;
B_an_n6 = ss_An.B;
C_an_n6 = ss_An.C;
D_an_n6 = ss_An.D;

W_an0_n6 = dcgain(W_an);

% annihilator neighbour 4
N = C_n4;
N_perp = null(N);
g = s+2;
W_an = zpk(minreal((N_perp*(s+1))/g));

ss_An = ss(W_an);

A_an_n4 = ss_An.A;
B_an_n4 = ss_An.B;
C_an_n4 = ss_An.C;
D_an_n4 = ss_An.D;

W_an0_n4 = dcgain(W_an);

% optimizer 

alpha = 50;
C_opt = eye(2);
x_opt0 = zeros(2,1);

% formation parameters

delta_n2 = [delta_x_t2(Idx2);delta_y_t2(Idx2)];
x0_n2 = x0_cc2(Idx2);
y0_n2 = y0_cc2(Idx2);
u_max_n2 = [u_max_x(idx_n2);u_max_y(idx_n2)];

delta_n3 = [delta_x_t2(Idx3);delta_y_t2(Idx3)];
x0_n3 = x0_cc2(Idx3);
y0_n3 = y0_cc2(Idx3);
u_max_n3 = [u_max_x(idx_n3);u_max_y(idx_n3)];

delta_n6 = [delta_x_t2(Idx6);delta_y_t2(Idx6)];
x0_n6 = x0_cc2(Idx6);
y0_n6 = y0_cc2(Idx6);
u_max_n6 = [u_max_x(idx_n6);u_max_y(idx_n6)];

delta_n4 = [delta_x_t2(Idx4);delta_y_t2(Idx4)];
x0_n4 = x0_cc2(Idx4);
y0_n4 = y0_cc2(Idx4);
u_max_n4 = [u_max_x(idx_n4);u_max_y(idx_n4)];

num = max(size(x0_n2));
sig_vec = ones(num,1);

% simulation

out_dead_t2 = sim('dec_consensus_dead');

t_dead_t2 = out_dead_t2.t;
dim_t2 = max(size(t_dead_t2));
x_dead_t2 = zeros(num,dim_t2);
y_dead_t2 = zeros(num,dim_t2);
u_sat_dead_t2 = zeros(2*num,dim_t2);

x_dead_t2(1,:) = out_dead_t2.y2(:,1,:);
y_dead_t2(1,:) = out_dead_t2.y2(:,2,:);

x_dead_t2(2,:) = out_dead_t2.y3(:,1,:);
y_dead_t2(2,:) = out_dead_t2.y3(:,2,:);

x_dead_t2(3,:) = out_dead_t2.y4(:,1,:);
y_dead_t2(3,:) = out_dead_t2.y4(:,2,:);

x_dead_t2(4,:) = out_dead_t2.y6(:,1,:);
y_dead_t2(4,:) = out_dead_t2.y6(:,2,:);

u_sat_dead_t2(1:2,:) = out_dead_t2.u_n2([2,5],:,:);
u_sat_dead_t2(3:4,:) = out_dead_t2.u_n3([2,5],:,:);
u_sat_dead_t2(5:6,:) = out_dead_t2.u_n4([2,5],:,:);
u_sat_dead_t2(7:8,:) = out_dead_t2.u_n6([2,5],:,:);

% figure()
% hold on; grid on;
% plot(x_dead_t2(:,1),y_dead_t2(:,1),'o');
% plot(x_dead_t2(:,end),y_dead_t2(:,end),'*');
% plot(x_dead_t2(1,:),y_dead_t2(1,:));
% plot(x_dead_t2(2,:),y_dead_t2(2,:));
% plot(x_dead_t2(3,:),y_dead_t2(3,:));
% plot(x_dead_t2(4,:),y_dead_t2(4,:));