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

A = zeros(n,n);

B = eye(n);

C = L;

% initial condition

x0 = zeros(n,1);

for i=1:n
    x0(i) = 8*i;
end

% saturation

u_max = [1,0.5,1,0.4,1.5,0.5,1]';

%% annihilator

% plant transfer function

s = tf('s');

Wp = C*(s*eye(n)-A)^(-1)*B;

% annihilator

N = L;

N_perp = null(L);

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

C_opt = 1;

x_opt0 = 0;

% threshold

th = 0.05;

%% time interval

dt = 0.01;
tf = 60;
tspan = 0:dt:tf;
dim = max(size(tspan));

f = 1;

%% nominal simulation

% array

x_nom = zeros(n,dim);

x_nom(:,1) = x0;

yc = zeros(n,dim);

% loop

delta = [2 1 0 -1 -2 -1  1]';

for i=1:dim-1   

    % current state
    xi = x_nom(:,i);
    yi = -L*(xi-delta) - sin(f*tspan(i));

    % controller
    yc(:,i) = yi;

    % euler integration
    x_nom(:,i+1) = x_nom(:,i) + dt*yi;
     

end

%% saturation simulation

% array

x_sat = zeros(n,dim);
x_sat(:,1) = x0;

yc_sat = zeros(n,dim);
u_sat = zeros(n,dim);

% loop

for i=1:dim-1

    % current state
    xi = x_sat(:,i);
    yi = -L*(xi-delta) - sin(f*tspan(i));

    % controller
    yc_sat(:,i) = yi;

    % saturation
    u_sat(:,i) = sat_fun(yc_sat(:,i),u_max);

    % euler integration
    x_sat(:,i+1) = xi + dt*u_sat(:,i);

end

%% deadzone optimization simulation

out_dead = sim('consensus_dead');

t_dead = out_dead.t;
x_dead = zeros(n,dim);
x_dead(:,:) = out_dead.y(:,:,:);
u_sat_dead = zeros(n,dim);
u_sat_dead(:,:) = out_dead.u(:,:,:);

%% transient optimization simulation

out_trans = sim('consensus_trans');

t_trans = out_trans.t;
x_trans = zeros(n,dim);
x_trans(:,:) = out_trans.y(:,:,:);
u_sat_trans = zeros(n,dim);
u_sat_trans(:,:) = out_trans.u(:,:,:);

%% nominal plot

figure()
hold on; grid on;
plot(tspan,x_nom)
title('consensus')

figure()
hold on; grid on;
plot(tspan,yc)
title('control')

figure()
hold on; grid on;
plot(tspan,L*(x_nom-delta))
title('error')


%% saturation plot

figure()
hold on; grid on;
plot(tspan,x_sat)
title('saturation consensus')

figure()
hold on; grid on;
plot(tspan,u_sat)
title('saturation control')

figure()
hold on; grid on;
plot(tspan,L*(x_sat-delta))
title('sat error')

%% deadzone plot

figure()
hold on; grid on;
plot(t_dead,x_dead)
title('dead consensus')

figure()
hold on; grid on;
plot(t_dead,u_sat_dead)
title('dead control')

figure()
hold on; grid on;
plot(t_dead,L*(x_dead-delta))
title('dead error')

%% transient plot

figure()
hold on; grid on;
plot(t_trans,x_trans)
title('trans consensus')

figure()
hold on; grid on;
plot(t_trans,u_sat_trans)
title('trans control')

figure()
hold on; grid on;
plot(t_trans,L*(x_trans-delta))
title('trans error')

%% functions

function u_sat = sat_fun(u,u_max)
    
    u_sat = u;
    for i=1:max(size(u))
        if (abs(u(i))>u_max(i))
            u_sat(i) = sign(u(i))*u_max(i);
        end
    end

end

function G = switch_gain(e,th)
    
    eps = th/2;

    G = 1;

    if ( norm(e) < th-eps )
    G = 0;
    end
    if ( norm(e) > th )
        G = 1;
    end
    if ( norm(e) <= th && norm(e) >= th-eps )
        G = -2*(norm(e)-(th-eps))^3+3*(norm(e)-(th-eps))^2;
    end

end






