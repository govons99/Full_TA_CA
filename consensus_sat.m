function [x_sat,y_sat,u_sat,tspan] = consensus_sat(L,x0,y0,dt,Tf,delta_x,delta_y,u_max)

% NOMINAL SIMULATION

tspan = 0:dt:Tf;
dim = max(size(tspan));

n = max(size(L));

% array

x_sat = zeros(n,dim);
y_sat = zeros(n,dim);

x_sat(:,1) = x0;
y_sat(:,1) = y0;

yc_sat = zeros(2*n,dim);
u_sat = zeros(2*n,dim);

% loop

f = 0.2;

for i=1:dim-1   

    % current state
    xi = x_sat(:,i);
    yi = y_sat(:,i);
    yci = [-L*(xi-delta_x) - cos(f*tspan(i));-L*(yi-delta_y) - sin(f*tspan(i))];

    % controller
    yc_sat(:,i) = yci;

    %saturation
    u_sat(:,i) = sat_fun(yc_sat(:,i),u_max);

    % euler integration
    x_sat(:,i+1) = x_sat(:,i) + dt*u_sat(1:n,i);
    y_sat(:,i+1) = y_sat(:,i) + dt*u_sat(n+1:end,i); 

end



