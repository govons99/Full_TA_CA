function [x_nom,y_nom,yc,tspan] = consensus(L,x0,y0,dt,Tf,delta_x,delta_y)

% NOMINAL SIMULATION

tspan = 0:dt:Tf;
dim = max(size(tspan));

n = max(size(L));

% array

x_nom = zeros(n,dim);
y_nom = zeros(n,dim);

x_nom(:,1) = x0;
y_nom(:,1) = y0;

yc = zeros(2*n,dim);

% loop

f = 0.2;

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



