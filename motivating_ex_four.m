clc
clear
close all

%% loading simulink block

addpath("./Simulink_motivating/")

%% time

% initial ref generator
dt1 = 0.01;
tspan = 0:dt1:1;

% robot simulation
dt = 0.05;

%% utilities

% number of robots

nr = 4;
nr2 = 2;
nr3 = 3;

% number of inputs

nu2 = 2*nr2;
nu3 = 2*nr3;

% number of states

nx2 = 2*nr2;
nx3 = 2*nr3;

% line parameter

a = 3;

% threshold for the tracking error

th = 0.05;

%% ref generation

vd = zeros(2*nr,max(size(tspan)));
pd = zeros(2*nr,max(size(tspan)));
pd(:,1) = [-2;4;4;4;4;0;-2;-2]; 

for i=1:max(size(tspan))-1
    
    p1 = pd(1:2,i);
    p2 = pd(3:4,i);
    p3 = pd(5:6,i);
    p4 = pd(7:8,i);
    
    sigma_l = [[a -1]*p1;[a -1]*p2;[a -1]*p3;[a -1]*p4];    
    Jl = blkdiag([a 1],[a 1],[a 1],[a 1]);    
    Jl_pinv = pinv(Jl);
    
    sigma_d = 0.5*[(p2-p1)'*(p2-p1);
                   (p3-p2)'*(p3-p2);
                   (p4-p3)'*(p4-p3)];
    Jd = 0.5*[-(p2-p1)',(p2-p1)',zeros(1,2),zeros(1,2);
              zeros(1,2),-(p3-p2)',(p3-p2)',zeros(1,2);
              zeros(1,2),zeros(1,2),-(p4-p3)',(p4-p3)'];
    Jd_pinv = pinv(Jd);
    
    sigma_tilde_l = zeros(nr,1) - sigma_l; 
    sigma_tilde_d = ones(nr-1,1) - sigma_d;
    
    v1 = Jl_pinv*10*sigma_tilde_l;
    v2 = Jd_pinv*6*sigma_tilde_d;
    
    vd(:,i+1) = v1 + (eye(2*nr)-Jl_pinv*Jl)*v2;
    
    pd(:,i+1) = pd(:,i)+vd(:,i+1)*dt1;
    
    
end

r_pd = pd(:,end);

t_change = 1000;
d = -4;

%% system

% static I/O linearization paramter

b = 0.2;

% output matrix: 2 agents
% x1 y1 x2 y2

C2 = [a -1 0 0;
      1 0 -1 0;
      0 0 a -1;
      0 1 0 -1
     ];

C_id2 = eye(nx2);

% output matrix: 3 agents
% x1 y1 x2 y2 x3 y3

C3 = [a -1 0 0 0 0;
      0 1 0 -1 0 0;
      0 0 a -1 0 0;
      0 0 1 0 -1 0;
      0 0 0 0 a -1;
      0 0 0 1 0 -1;
      ];

C_id3 = eye(nx3);
 
% plant

s = tf('s');

A2 = zeros(nx2);
B2 = eye(nu2);

A3 = zeros(nx3);
B3 = eye(nu3);

Wp2 = C2*inv(s*eye(nx2)-A2)*B2;
Wp3 = C3*inv(s*eye(nx3)-A3)*B3;

Wp_id2 = C_id2*inv(s*eye(nx2)-A2)*B2;
Wp_id3 = C_id3*inv(s*eye(nx3)-A3)*B3;

%% controller redundancy

Wc1 = pidtune(Wp2(1,1),'PID');
Wc2 = pidtune(Wp2(2,1),'PID');
Wc3 = pidtune(Wp2(3,3),'PID');
Wc4 = pidtune(Wp2(4,4),'PID');

Wc_mrs2 = blkdiag(Wc1,Wc2,Wc3,Wc4)*2;

ssC_mrs2 = ss(Wc_mrs2);

Ac_mrs2 = ssC_mrs2.A;
Bc_mrs2 = ssC_mrs2.B;
Cc_mrs2 = ssC_mrs2.C;
Dc_mrs2 = ssC_mrs2.D;

Wc1 = pidtune(Wp3(1,1),'PID');
Wc2 = pidtune(Wp3(2,2),'PID');
Wc3 = pidtune(Wp3(3,3),'PID');
Wc4 = pidtune(Wp3(4,3),'PID');
Wc5 = pidtune(Wp3(5,5),'PID');
Wc6 = pidtune(Wp3(6,6),'PID');

Wc_mrs3 = blkdiag(Wc1,Wc2,Wc3,Wc4,Wc5,Wc6)*2;

ssC_mrs3 = ss(Wc_mrs3);

Ac_mrs3 = ssC_mrs3.A;
Bc_mrs3 = ssC_mrs3.B;
Cc_mrs3 = ssC_mrs3.C;
Dc_mrs3 = ssC_mrs3.D;

Kp = 0.85;
Ki = 1;

%% allocator: annihilator + optimizer

% annihilator

N2 = C2;
N3 = C3;

N_perp2 = null(C2);
N_perp3 = null(C3);

g = s+2;

W_an2 = zpk(minreal((N_perp2*(s+1))/g));

ss_An2 = ss(W_an2);

A_an2 = ss_An2.A;
B_an2 = ss_An2.B;
C_an2 = ss_An2.C;
D_an2 = ss_An2.D;

W_an02 = dcgain(W_an2);

xan02 = zeros(max(size(A_an2)),1);

W_an3 = zpk(minreal((N_perp3*(s+1))/g));

ss_An3 = ss(W_an3);

A_an3 = ss_An3.A;
B_an3 = ss_An3.B;
C_an3 = ss_An3.C;
D_an3 = ss_An3.D;

W_an03 = dcgain(W_an3);

xan03 = zeros(max(size(A_an3)),1);

% optimizer

alpha = 25;

C_opt = 1;

xopt0 = zeros(max(size(C_opt)),1);

%% saturation limits

% u_max1 = [0.75;1.25];
% u_max2 = [0.70;1.15];
% u_max3 = [0.5;1.15];
% u_max4 = [0.5;1.25];
% u_max5 = [0.5;1.25];
% u_max6 = [0.5;1.25];

% u_max1 = [0.75;1.15];
% u_max2 = [0.70;1.15];
% u_max3 = [0.5;1.15];
% u_max4 = [0.5;1.15];
% u_max5 = [0.5;1.25];
% u_max6 = [0.5;1.25];

u_max1 = [0.3;0.5];
u_max2 = [0.15;0.3];
u_max3 = [0.15;0.3];
u_max4 = [0.15;0.3];

u_max_mrs = [u_max1;u_max2;u_max3;u_max4];

u_max_n1 = [u_max1;u_max2];
u_max_n2 = [u_max1;u_max2;u_max3];
u_max_n3 = [u_max2;u_max3;u_max4];
u_max_n4 = [u_max3;u_max4];

%% robot simulation

%[-2;4;4;4;4;0;-2;-2;-4;-4;1;-4]

x10 = -2; y10 = 4;
x20 = 4;  y20 = 4;
x30 = 4;  y30 = 0;
x40 = -2; y40 = -2;

p01 = [x10;y10];
p02 = [x20;y20];
p03 = [x30;y30];
p04 = [x40;y40];

% p01 = r_pd(1:2);
% p02 = r_pd(3:4);
% p03 = r_pd(5:6);
% p04 = r_pd(7:8);
% p05 = r_pd(9:10);
% p06 = r_pd(11:12);

th02 = atan(a)*ones(nr2,1);
phi02 = atan(a)*ones(nr2,1);

th03 = atan(a)*ones(nr3,1);
phi03 = atan(a)*ones(nr3,1);

xc02 = zeros(max(size(Ac_mrs2)),1);
xc03 = zeros(max(size(Ac_mrs3)),1);

%% nominal condition

Tf = 150;

% first neighbour N1

p_in1 = [p01;p02];
r_pd1 = r_pd(1:4);

% second neighbout N2

p_in2 = [p01;p02;p03];
r_pd2 = r_pd(1:6);

% third neighbour N3

p_in3 = [p02;p03;p04];
r_pd3 = r_pd(3:8);

% fourth neighbour N4

p_in4 = [p03;p04];
r_pd4 = r_pd(5:8);

out = sim('decentralized');

time = out.t;

r1 = out.y1(:,1:2);
theta1 = out.theta1(:,1);
u1 = out.u1(:,1:2);

r2 = out.y2(:,1:2);
theta2 = out.theta2(:,2);
u2 = out.u2(:,3:4);

r3 = out.y3(:,1:2);
theta3 = out.theta3(:,2);
u3 = out.u3(:,3:4);

r4 = out.y4(:,1:2);
theta4 = out.theta4(:,2);
u4 = out.u4(:,3:4);

%% saturating condition

out_sat = sim('decentralized_sat');

time_sat = out_sat.t;

r1_sat = out_sat.y1(:,1:2);
theta1_sat = out_sat.theta1(:,1);
u1_sat = out_sat.u1(:,1:2);

r2_sat = out_sat.y2(:,1:2);
theta2_sat = out_sat.theta2(:,2);
u2_sat = out_sat.u2(:,3:4);

r3_sat = out_sat.y3(:,1:2);
theta3_sat = out_sat.theta3(:,2);
u3_sat = out_sat.u3(:,3:4);

r4_sat = out_sat.y4(:,1:2);
theta4_sat = out_sat.theta4(:,2);
u4_sat = out_sat.u4(:,3:4);

%% minimization of the deadzone at steady-state

out_dead = sim('decentralized_dead');

time_dead = out_dead.t;

r1_dead = out_dead.y1(:,1:2);
theta1_dead = out_dead.theta1(:,1);
u1_dead = out_dead.u1(:,1:2);

r2_dead = out_dead.y2(:,1:2);
theta2_dead = out_dead.theta2(:,2);
u2_dead = out_dead.u2(:,3:4);

r3_dead = out_dead.y3(:,1:2);
theta3_dead = out_dead.theta3(:,2);
u3_dead = out_dead.u3(:,3:4);

r4_dead = out_dead.y4(:,1:2);
theta4_dead = out_dead.theta4(:,2);
u4_dead = out_dead.u4(:,3:4);

% %% minimization of the deadzone along the transient
% 
% % t1 = 45;
% % t2 = 85;
% % t3 = 115;
% % Tf = 150;
% 
% out_trans = sim('decentralized_trans');
% 
% time_trans = out_trans.t;
% 
% r1_trans = out_trans.y1(:,1:2);
% theta1_trans = out_trans.theta1(:,1);
% u1_trans = out_trans.u1(:,1:2);
% 
% r2_trans = out_trans.y2(:,1:2);
% theta2_trans = out_trans.theta2(:,2);
% u2_trans = out_trans.u2(:,3:4);
% 
% r3_trans = out_trans.y3(:,1:2);
% theta3_trans = out_trans.theta3(:,2);
% u3_trans = out_trans.u3(:,3:4);
% 
% r4_trans = out_trans.y4(:,1:2);
% theta4_trans = out_trans.theta4(:,2);
% u4_trans = out_trans.u4(:,3:4);
% 
% r5_trans = out_trans.y5(:,1:2);
% theta5_trans = out_trans.theta5(:,2);
% u5_trans = out_trans.u5(:,3:4);
% 
% r6_trans = out_trans.y6(:,1:2);
% theta6_trans = out_trans.theta6(:,2);
% u6_trans = out_trans.u6(:,3:4);

%% workspace

d = datetime("today");
%save(['decentralized_',d]);

%% nominal plots

figure()
hold on; grid on;
p1 = plot(r1(:,1),r1(:,2));
p2 = plot(r2(:,1),r2(:,2));
p3 = plot(r3(:,1),r3(:,2));
p4 = plot(r4(:,1),r4(:,2));
plot(r1(end,1),r1(end,2),'o','linewidth',2);
plot(r2(end,1),r2(end,2),'o','linewidth',2);
plot(r3(end,1),r3(end,2),'o','linewidth',2);
plot(r4(end,1),r4(end,2),'o','linewidth',2);
line([-3,2],[-9,6],'linewidth',1.5)
legend([p1,p2,p3,p4],'r1','r2','r3','r4','location','nw');
%saveas(gcf,'./dec/platoon_nominal_traj.eps','epsc')

figure()

tiledlayout(1,4)

nexttile
hold on; grid on;
plot(time,u1,'linewidth',2);
title('robot 1');
xlim([0,time(end)]);
lg = legend('$u_1$','$u_2$','numcolumns',2,'interpreter','latex');
lg.Layout.Tile = 'north';
ylabel("$u [m/s^2]$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time,u2,'linewidth',2);
title('robot 2');
xlim([0,time(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time,u3,'linewidth',2);
title('robot 3');
xlim([0,time(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time,u4,'linewidth',2);
title('robot 4');
xlim([0,time(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

% nexttile
% hold on; grid on;
% plot(time,u4,'linewidth',2);
% title('robot 4');
% xlim([0,time(end)]);
% ylabel("$u [m/s^2]$",'interpreter','latex')
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;
% 
% nexttile
% hold on; grid on;
% plot(time,u5,'linewidth',2);
% title('robot 5');
% xlim([0,time(end)]);
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;
% 
% nexttile
% hold on; grid on;
% plot(time,u3,'linewidth',2);
% title('robot 6');
% xlim([0,time(end)]);
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;

%saveas(gcf,'./dec/platoon_nominal_control.eps','epsc');

%% saturating plots

figure()
hold on; grid on;
p1 = plot(r1_sat(:,1),r1_sat(:,2),'linewidth',2);
p2 = plot(r2_sat(:,1),r2_sat(:,2),'linewidth',2);
p3 = plot(r3_sat(:,1),r3_sat(:,2),'linewidth',2);
p4 = plot(r4_sat(:,1),r4_sat(:,2),'linewidth',2);
plot(r1_sat(end,1),r1_sat(end,2),'o','linewidth',2);
plot(r2_sat(end,1),r2_sat(end,2),'o','linewidth',2);
plot(r3_sat(end,1),r3_sat(end,2),'o','linewidth',2);
plot(r4_sat(end,1),r4_sat(end,2),'o','linewidth',2);
% line([-4,2],[-12,6],'linewidth',1.5)
ylabel("$y$",'interpreter','latex')
xlabel("$x$",'interpreter','latex')
legend([p1,p2,p3,p4],'$r_1$','$r_2$','$r_3$','$r_4$','location','nw','interpreter','latex');
set(gca,'FontSize',16)
saveas(gcf,'./dec/platoon_sat_traj.eps','epsc')

figure()

tiledlayout(1,4)

nexttile
hold on; grid on;
plot(time_sat,u1_sat,'linewidth',2);
title('robot 1');
xlim([0,time_sat(end)]);
lg = legend('$u_1$','$u_2$','numcolumns',2,'interpreter','latex');
lg.Layout.Tile = 'north';
ylabel("$u [m/s^2]$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_sat,u2_sat,'linewidth',2);
title('robot 2');
xlim([0,time_sat(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_sat,u3_sat,'linewidth',2);
title('robot 3');
xlim([0,time_sat(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_sat,u4_sat,'linewidth',2);
title('robot 4');
xlim([0,time_sat(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;


% nexttile
% hold on; grid on;
% plot(time_sat,u4_sat,'linewidth',2);
% title('robot 4');
% xlim([0,time_sat(end)]);
% ylabel("$u [m/s^2]$",'interpreter','latex')
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;
% 
% nexttile
% hold on; grid on;
% plot(time_sat,u5_sat,'linewidth',2);
% title('robot 5');
% xlim([0,time_sat(end)]);
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;
% 
% nexttile
% hold on; grid on;
% plot(time_sat,u6_sat,'linewidth',2);
% title('robot 6');
% xlim([0,time_sat(end)]);
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;

%saveas(gcf,'./dec/platoon_sat_control.eps','epsc');

%% minimization of the deadzone at steady-state

% trajectory
figure()
hold on; grid on;
p1 = plot(r1_dead(:,1),r1_dead(:,2),'linewidth',2);
p2 = plot(r2_dead(:,1),r2_dead(:,2),'linewidth',2);
p3 = plot(r3_dead(:,1),r3_dead(:,2),'linewidth',2);
p4 = plot(r4_dead(:,1),r4_dead(:,2),'linewidth',2);
plot(r1_dead(end,1),r1_dead(end,2),'o','linewidth',2);
plot(r2_dead(end,1),r2_dead(end,2),'o','linewidth',2);
plot(r3_dead(end,1),r3_dead(end,2),'o','linewidth',2);
plot(r4_dead(end,1),r4_dead(end,2),'o','linewidth',2);
%line([-4,6],[-12,18],'linewidth',1.5)
xlim([-3,6]);
ylim([-4,12]);
ylabel("$y$",'interpreter','latex')
xlabel("$x$",'interpreter','latex')
%title('all dead')
legend([p1,p2,p3,p4],'$r_1$','$r_2$','$r_3$','$r_4$','location','nw','interpreter','latex');
set(gca,'FontSize',16)
saveas(gcf,'./dec/platoon_dead_traj.eps','epsc')

% total control
figure()

tiledlayout(1,4)

nexttile
hold on; grid on;
plot(time_dead,u1_dead,'linewidth',2);
title('robot 1');
xlim([0,time_dead(end)]);
lg = legend('$u_1$','$u_2$','numcolumns',2,'interpreter','latex');
lg.Layout.Tile = 'north';
ylabel("$u [m/s^2]$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_dead,u2_dead,'linewidth',2);
title('robot 2');
xlim([0,time_dead(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_dead,u3_dead,'linewidth',2);
title('robot 3');
xlim([0,time_dead(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_dead,u4_dead,'linewidth',2);
title('robot 4');
xlim([0,time_dead(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

% nexttile
% hold on; grid on;
% plot(time_dead,u4_dead,'linewidth',2);
% title('robot 4');
% xlim([0,time_dead(end)]);
% ylabel("$u [m/s^2]$",'interpreter','latex')
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;
% 
% nexttile
% hold on; grid on;
% plot(time_dead,u5_dead,'linewidth',2);
% title('robot 5');
% xlim([0,time_dead(end)]);
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;
% 
% nexttile
% hold on; grid on;
% plot(time_dead,u6_dead,'linewidth',2);
% title('robot 6');
% xlim([0,time_dead(end)]);
% xlabel("$Time [s]$",'interpreter','latex')
% hold off;

%saveas(gcf,'./dec/platoon_dead_control.eps','epsc');

% nominal controller
% figure()
% 
% subplot(221)
% hold on; grid on;
% plot(time_dead,u1_dead,'linewidth',2);
% title('control robot 1');
% xlim([0,Tf]);
% 
% subplot(222)
% hold on; grid on;
% plot(time_dead,u2_dead,'linewidth',2);
% title('control robot 2');
% xlim([0,Tf]);
% 
% subplot(223)
% hold on; grid on;
% plot(time_dead,u3_dead,'linewidth',2);
% title('control robot 3');
% xlim([0,Tf]);
% 
% subplot(224)
% hold on; grid on;
% plot(time_dead,u4_dead,'linewidth',2);
% title('control robot 4');
% xlim([0,Tf]);
% 
% sgtitle('all dead nominal control')

%% minimization of the deadzone along the transient

figure()
hold on; grid on;
p1 = plot(r1_trans(:,1),r1_trans(:,2));
p2 = plot(r2_trans(:,1),r2_trans(:,2));
p3 = plot(r3_trans(:,1),r3_trans(:,2));
p4 = plot(r4_trans(:,1),r4_trans(:,2));
p5 = plot(r5_trans(:,1),r5_trans(:,2));
p6 = plot(r6_trans(:,1),r6_trans(:,2));
plot(r1_trans(end,1),r1_trans(end,2),'o','linewidth',2);
plot(r2_trans(end,1),r2_trans(end,2),'o','linewidth',2);
plot(r3_trans(end,1),r3_trans(end,2),'o','linewidth',2);
plot(r4_trans(end,1),r4_trans(end,2),'o','linewidth',2);
plot(r5_trans(end,1),r5_trans(end,2),'o','linewidth',2);
plot(r6_trans(end,1),r6_trans(end,2),'o','linewidth',2);
line([-3,4],[-9,12],'linewidth',1.5)
title('all trans')
legend([p1,p2,p3,p4,p5,p6],'r1','r2','r3','r4','r5','r6','location','nw');
%%saveas(gcf,'./dec/platoon_trans_traj.eps','epsc')

% total control
figure()

tiledlayout(2,3)

nexttile
hold on; grid on;
plot(time_trans,u1_trans,'linewidth',2);
title('robot 1');
xlim([0,time_trans(end)]);
lg = legend('$u_1$','$u_2$','numcolumns',2,'interpreter','latex');
lg.Layout.Tile = 'north';
ylabel("$u [m/s^2]$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_trans,u2_trans,'linewidth',2);
title('robot 2');
xlim([0,time_trans(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_trans,u3_trans,'linewidth',2);
title('robot 3');
xlim([0,time_trans(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_trans,u4_trans,'linewidth',2);
title('robot 4');
xlim([0,time_trans(end)]);
ylabel("$u [m/s^2]$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_trans,u5_trans,'linewidth',2);
title('robot 5');
xlim([0,time_trans(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_trans,u6_trans,'linewidth',2);
title('robot 6');
xlim([0,time_trans(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

%%saveas(gcf,'./dec/platoon_trans_control.eps','epsc');

% nominal control
% figure()
% 
% subplot(221)
% hold on; grid on;
% plot(time_trans,u1_trans,'linewidth',2);
% title('control robot 1');
% xlim([0,Tf]);
% 
% subplot(222)
% hold on; grid on;
% plot(time_trans,u2_trans,'linewidth',2);
% title('control robot 2');
% xlim([0,Tf]);
% 
% subplot(223)
% hold on; grid on;
% plot(time_trans,u3_trans,'linewidth',2);
% title('control robot 3');
% xlim([0,Tf]);
% 
% subplot(224)
% hold on; grid on;
% plot(time_trans,u4_trans,'linewidth',2);
% title('control robot 4');
% xlim([0,Tf]);
% 
% sgtitle('all trans nominal control')

%% 

figure()

tiledlayout(2,4)

nexttile
hold on; grid on;
plot(time_sat,u1_sat,'linewidth',2);
title('robot 1');
xlim([0,time_sat(end)]);
lg = legend('$u_1$','$u_2$','numcolumns',2,'interpreter','latex');
lg.Layout.Tile = 'north';
ylabel("$u [m/s^2]$",'interpreter','latex')
%xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_sat,u2_sat,'linewidth',2);
title('robot 2');
xlim([0,time_sat(end)]);
%xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_sat,u3_sat,'linewidth',2);
title('robot 3');
xlim([0,time_sat(end)]);
%xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_sat,u4_sat,'linewidth',2);
title('robot 4');
xlim([0,time_sat(end)]);
%xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_dead,u1_dead,'linewidth',2);
%title('robot 1');
xlim([0,time_dead(end)]);
% lg = legend('$u_1$','$u_2$','numcolumns',2,'interpreter','latex');
% lg.Layout.Tile = 'north';
ylabel("$u [m/s^2]$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_dead,u2_dead,'linewidth',2);
%title('robot 2');
xlim([0,time_dead(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_dead,u3_dead,'linewidth',2);
%title('robot 3');
xlim([0,time_dead(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

nexttile
hold on; grid on;
plot(time_dead,u4_dead,'linewidth',2);
%title('robot 3');
xlim([0,time_dead(end)]);
xlabel("$Time [s]$",'interpreter','latex')
hold off;

%saveas(gcf,'./dec/comp_control.eps','epsc');
    
    