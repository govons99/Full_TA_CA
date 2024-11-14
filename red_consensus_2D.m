clc
clear
close all

%% graph matrices

ad = [0, 0, 0, 0, 1, 0, 0;
      0, 0, 1, 0, 1, 0, 0;
      0, 1, 0, 1, 0, 0, 1;
      0, 0, 1, 0, 0, 1, 0;
      1, 1, 0, 0, 0, 0, 0;
      0, 0, 0, 1, 0, 0, 0;
      0, 0, 1, 0, 0, 0, 0];

D = diag([1,2,3,2,2,1,1]);

L = D - ad;

%% single integrator dynamics

% number of robots

n = max(size(ad));

% dynamical matrices

A = zeros(2*n,2*n);

B = eye(2*n);

C = blkdiag(L,L);

% initial condition

x0 = zeros(n,1);
y0 = zeros(n,1);

for i=1:n
    x0(i) = 8*i;
    y0(i) = 8*i;
end

% saturation

u_max_x = [1,0.5,1,0.4,1.5,0.5,1]';
u_max_y = [1,0.5,1,0.4,1.5,0.5,1]';
u_max = [u_max_x;u_max_y];

%% annihilator

% plant transfer function

s = tf('s');

Wp = C*(s*eye(2*n)-A)^(-1)*B;

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

%% time interval

dt = 0.01;
tf = 100;
tspan = 0:dt:tf;
dim = max(size(tspan));

f = 0.2;

%% nominal simulation

% array

x_nom = zeros(n,dim);
y_nom = zeros(n,dim);

x_nom(:,1) = x0;
y_nom(:,1) = y0;

yc = zeros(2*n,dim);

% loop

delta_x = [2 0 1 -1 -2 -1  1]';
delta_y = [0 0 2  2  0 -2 -2]';

% delta_x = zeros(n,1);
% delta_y = zeros(n,1);

delta =  [delta_x;delta_y];

for i=1:dim-1   

    % current state
    xi = x_nom(:,i);
    yi = y_nom(:,i);
    yci = [-L*(xi-delta_x) - cos(f*tspan(i));-L*(yi-delta_y) - sin(f*tspan(i))];

    % controller
    yc(:,i) = yci;

    % euler integration
    x_nom(:,i+1) = x_nom(:,i) + dt*yci(1:n);
    y_nom(:,i+1) = y_nom(:,i) + dt*yci(n+1:end); 

end

%% saturation simulation

% array

x_sat = zeros(n,dim);
y_sat = zeros(n,dim);
x_sat(:,1) = x0;
y_sat(:,1) = y0;

yc_sat = zeros(2*n,dim);
u_sat = zeros(2*n,dim);

% loop

for i=1:dim-1

    % current state
    xi = x_sat(:,i);
    yi = y_sat(:,i);
    yci = [-L*(xi-delta_x) - cos(f*tspan(i));-L*(yi-delta_y) - sin(f*tspan(i))];

    % controller
    yc_sat(:,i) = yci;

    % saturation
    u_sat(:,i) = sat_fun(yc_sat(:,i),u_max);

    % euler integration
    x_sat(:,i+1) = xi + dt*u_sat(1:n,i);
    y_sat(:,i+1) = yi + dt*u_sat(n+1:end,i);

end

%% deadzone optimization simulation

out_dead = sim('consensus_dead');

t_dead = out_dead.t;
x_dead = zeros(n,dim);
y_dead = zeros(n,dim);
x_dead(:,:) = out_dead.y(1:n,:,:);
y_dead(:,:) = out_dead.y(n+1:end,:,:);
u_sat_dead = zeros(2*n,dim);
u_sat_dead(:,:) = out_dead.u(:,:,:);

%% transient optimization simulation

out_trans = sim('consensus_trans');

t_trans = out_trans.t;
x_trans = zeros(n,dim);
y_trans = zeros(n,dim);
x_trans(:,:) = out_trans.y(1:n,:,:);
y_trans(:,:) = out_trans.y(n+1:end,:,:);
u_sat_trans = zeros(2*n,dim);
u_sat_trans(:,:) = out_trans.u(:,:,:);

%% nominal plot

figure()
hold on; grid on;
plot(x_nom(:,1),y_nom(:,1),'o');
plot(x_nom(:,end),y_nom(:,end),'*');
plot(x_nom(1,:),y_nom(1,:));
plot(x_nom(2,:),y_nom(2,:));
plot(x_nom(3,:),y_nom(3,:));
plot(x_nom(4,:),y_nom(4,:));
plot(x_nom(5,:),y_nom(5,:));
plot(x_nom(6,:),y_nom(6,:));
plot(x_nom(7,:),y_nom(7,:));

figure()

subplot(211)
hold on; grid on;
plot(tspan,x_nom)
title('x-consensus')

subplot(212)
hold on; grid on;
plot(tspan,y_nom)
title('y-consensus')

figure()
hold on; grid on;
plot(tspan,yc)
title('control')

figure()

subplot(211)
hold on; grid on;
plot(tspan,L*(x_nom-delta_x))
title('x-error')

subplot(212)
hold on; grid on;
plot(tspan,L*(y_nom-delta_y))
title('y-error')


%% saturation plot

figure()
hold on; grid on;
plot(x_sat(:,1),y_sat(:,1),'o');
plot(x_sat(:,end),y_sat(:,end),'*');
plot(x_sat(1,:),y_sat(1,:));
plot(x_sat(2,:),y_sat(2,:));
plot(x_sat(3,:),y_sat(3,:));
plot(x_sat(4,:),y_sat(4,:));
plot(x_sat(5,:),y_sat(5,:));
plot(x_sat(6,:),y_sat(6,:));
plot(x_sat(7,:),y_sat(7,:));

figure()

subplot(211)
hold on; grid on;
plot(tspan,x_sat)
title('saturation x-consensus')

subplot(212)
hold on; grid on;
plot(tspan,y_sat)
title('saturation y-consensus')

figure()

subplot(211)
hold on; grid on;
plot(tspan,u_sat(1:n,:))
title('saturation x-control')

subplot(212)
hold on; grid on;
plot(tspan,u_sat(n+1:end,:))
title('saturation y-control')

figure()

subplot(211)
hold on; grid on;
plot(tspan,L*(x_sat-delta_x))
title('sat x-error')

subplot(212)
hold on; grid on;
plot(tspan,L*(y_sat-delta_y))
title('sat y-error')

%% deadzone plot

figure()
hold on; grid on;
plot(x_dead(:,1),y_dead(:,1),'o');
plot(x_dead(:,end),y_dead(:,end),'*');
plot(x_dead(1,:),y_dead(1,:));
plot(x_dead(2,:),y_dead(2,:));
plot(x_dead(3,:),y_dead(3,:));
plot(x_dead(4,:),y_dead(4,:));
plot(x_dead(5,:),y_dead(5,:));
plot(x_dead(6,:),y_dead(6,:));
plot(x_dead(7,:),y_dead(7,:));

figure()

subplot(211)
hold on; grid on;
plot(t_dead,x_dead)
title('dead x-consensus')

subplot(212)
hold on; grid on;
plot(t_dead,y_dead)
title('dead y-consensus')

figure()

subplot(211)
hold on; grid on;
plot(t_dead,u_sat_dead(1:n,:))
title('dead x-control')

subplot(212)
hold on; grid on;
plot(t_dead,u_sat_dead(n+1:end,:))
title('dead y-control')

figure()

subplot(211)
hold on; grid on;
plot(t_dead,L*(x_dead-delta_x))
title('dead x-error')

subplot(212)
hold on; grid on;
plot(t_dead,L*(y_dead-delta_y))
title('dead y-error')

%% transient plot

figure()
hold on; grid on;
plot(x_trans(:,1),y_trans(:,1),'o');
plot(x_trans(:,end),y_trans(:,end),'*');
plot(x_trans(1,:),y_trans(1,:));
plot(x_trans(2,:),y_trans(2,:));
plot(x_trans(3,:),y_trans(3,:));
plot(x_trans(4,:),y_trans(4,:));
plot(x_trans(5,:),y_trans(5,:));
plot(x_trans(6,:),y_trans(6,:));
plot(x_trans(7,:),y_trans(7,:));

figure()

subplot(211)
hold on; grid on;
plot(t_trans,x_trans)
title('trans x-consensus')

subplot(212)
hold on; grid on;
plot(t_trans,y_trans)
title('trans y-consensus')

figure()

subplot(211)
hold on; grid on;
plot(t_trans,u_sat_trans(1:n,:))
title('trans x-control')

subplot(212)
hold on; grid on;
plot(t_trans,u_sat_trans(n+1:end,:))
title('trans control')

figure()

subplot(211)
hold on; grid on;
plot(t_trans,L*(x_trans-delta_x))
title('trans x-error')

subplot(212)
hold on; grid on;
plot(t_trans,L*(x_trans-delta_y))
title('trans y-error')


