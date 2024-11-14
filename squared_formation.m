% execution of the squared formation

delta_x_t2 = [0, 0, 2, 2]';
delta_y_t2 = [0, 2, 0, 2]';

ad_t2 = ad(cc{idx2},cc{idx2});

D_t2 = diag(sum(ad_t2,1));

Lt2 = D_t2 - ad_t2;

u_max = [u_max_x(cc{idx2});u_max_y(cc{idx2})];

[x_nom_t2,y_nom_t2,yc_t2,tspan_t2] = consensus(Lt2,x0_cc2,y0_cc2,dt,Tf,delta_x_t2,delta_y_t2);

% figure()
% hold on; grid on;
% plot(x_nom_t2(:,1),y_nom_t2(:,1),'o');
% plot(x_nom_t2(:,end),y_nom_t2(:,end),'*');
% plot(x_nom_t2(1,:),y_nom_t2(1,:));
% plot(x_nom_t2(2,:),y_nom_t2(2,:));
% plot(x_nom_t2(3,:),y_nom_t2(3,:));
% plot(x_nom_t2(4,:),y_nom_t2(4,:));

[x_sat_t2,y_sat_t2,yc_sat_t2,tspan_sat_t2] = consensus_sat(Lt2,x0_cc2,y0_cc2,dt,Tf,delta_x_t2,delta_y_t2,u_max);

% figure()
% hold on; grid on;
% plot(x_sat_t2(:,1),y_sat_t2(:,1),'o');
% plot(x_sat_t2(:,end),y_sat_t2(:,end),'*');
% plot(x_sat_t2(1,:),y_sat_t2(1,:));
% plot(x_sat_t2(2,:),y_sat_t2(2,:));
% plot(x_sat_t2(3,:),y_sat_t2(3,:));
% plot(x_sat_t2(4,:),y_sat_t2(4,:));

% AGENT DYNAMICS

% number of robots
n_t2 = max(size(Lt2));

% dynamical matrices
A = zeros(2*n_t2,2*n_t2);
B = eye(2*n_t2);
C = blkdiag(Lt2,Lt2);

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

delta = [delta_x_t2;delta_y_t2];
x0 = x0_cc2;
y0 = y0_cc2;
num = max(size(x0));
sig_vec = ones(num,1);

out_dead_t2 = sim('consensus_dead');

t_dead_t2 = out_dead_t2.t;
dim_t2 = max(size(t_dead_t2));
x_dead_t2 = zeros(num,dim_t2);
y_dead_t2 = zeros(num,dim_t2);
x_dead_t2(:,:) = out_dead_t2.y(1:num,:,:);
y_dead_t2(:,:) = out_dead_t2.y(num+1:end,:,:);
u_sat_dead_t2 = zeros(2*num,dim_t2);
u_sat_dead_t2(:,:) = out_dead_t2.u(:,:,:);

% figure()
% hold on; grid on;
% plot(x_dead_t2(:,1),y_dead_t2(:,1),'o');
% plot(x_dead_t2(:,end),y_dead_t2(:,end),'*');
% plot(x_dead_t2(1,:),y_dead_t2(1,:));
% plot(x_dead_t2(2,:),y_dead_t2(2,:));
% plot(x_dead_t2(3,:),y_dead_t2(3,:));
% plot(x_dead_t2(4,:),y_dead_t2(4,:));