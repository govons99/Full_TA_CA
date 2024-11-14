% execution of the rendez-vous

delta_x_t3 = zeros(max(size(x0_cc3)),1);
delta_y_t3 = zeros(max(size(x0_cc3)),1);

ad_t3 = ad(cc{idx3},cc{idx3});

D_t3 = diag(sum(ad_t3,1));

Lt3 = D_t3 - ad_t3;

u_max = [u_max_x(cc{idx3});u_max_y(cc{idx3})];

[x_nom_t3,y_nom_t3,yc_t3,tspan_t3] = consensus(Lt3,x0_cc3,y0_cc3,dt,Tf,delta_x_t3,delta_y_t3);

% figure()
% hold on; grid on;
% plot(x_nom_t3(:,1),y_nom_t3(:,1),'o');
% plot(x_nom_t3(:,end),y_nom_t3(:,end),'*');
% plot(x_nom_t3(1,:),y_nom_t3(1,:));
% plot(x_nom_t3(2,:),y_nom_t3(2,:));
% plot(x_nom_t3(3,:),y_nom_t3(3,:));

[x_sat_t3,y_sat_t3,yc_sat_t3,tspan_sat_t3] = consensus_sat(Lt3,x0_cc3,y0_cc3,dt,Tf,delta_x_t3,delta_y_t3,u_max);

% figure()
% hold on; grid on;
% plot(x_sat_t3(:,1),y_sat_t3(:,1),'o');
% plot(x_sat_t3(:,end),y_sat_t3(:,end),'*');
% plot(x_sat_t3(1,:),y_sat_t3(1,:));
% plot(x_sat_t3(2,:),y_sat_t3(2,:));
% plot(x_sat_t3(3,:),y_sat_t3(3,:));

% AGENT DYNAMICS

% number of robots
n_t3 = max(size(Lt3));

% dynamical matrices
A = zeros(2*n_t3,2*n_t3);
B = eye(2*n_t3);
C = blkdiag(Lt3,Lt3);

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

f = 0.2;
delta = [delta_x_t3;delta_y_t3];
x0 = x0_cc3;
y0 = y0_cc3;
num = max(size(x0));
sig_vec = ones(num,1);

out_dead_t3 = sim('consensus_dead');

t_dead_t3 = out_dead_t3.t;
dim_t3 = max(size(t_dead_t3));
x_dead_t3 = zeros(num,dim_t3);
y_dead_t3 = zeros(num,dim_t3);
x_dead_t3(:,:) = out_dead_t3.y(1:num,:,:);
y_dead_t3(:,:) = out_dead_t3.y(num+1:end,:,:);
u_sat_dead_t3 = zeros(2*num,dim_t3);
u_sat_dead_t3(:,:) = out_dead_t3.u(:,:,:);

% figure()
% hold on; grid on;
% plot(x_dead_t3(:,1),y_dead_t3(:,1),'o');
% plot(x_dead_t3(:,end),y_dead_t3(:,end),'*');
% plot(x_dead_t3(1,:),y_dead_t3(1,:));
% plot(x_dead_t3(2,:),y_dead_t3(2,:));
% plot(x_dead_t3(3,:),y_dead_t3(3,:));
