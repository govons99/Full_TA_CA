function [z,d] = Task_All_new(nr,Pi,alpha,nt,R,delta,F,T,n_min,n_max,p,c,l,M)

ops = sdpsettings('verbose',0);

cost1 = 0;
cost2 = 0;
cost3 = 0;

% constraint vector

cons_vec = [];

% robot loop

for i=1:nr

    % bad allocation cost term
    temp = norm(Pi{i}*alpha(:,i))^2;
    cost1 = cost1+temp;

    % constraint b
    Bi = ones(1,nt)*alpha(:,i)<=1;
    cons_vec = vertcat(cons_vec,Bi);

    % constraint f
    sumi = sum(alpha(:,i));
    Fi = ( delta(i) <= sumi ) <= nt + (1-nt)*delta(i);
    cons_vec = vertcat(cons_vec,Fi);

    %constraint g
    for m=1:nt
        sumi = sum(alpha(:,i))-alpha(m,i);
        G_count = delta(i) >= alpha(m,i) - sumi;
        cons_vec = vertcat(cons_vec,G_count);
    end

end

% task loop

for m=1:nt

    % redundancy cost term
    temp = norm(R{m}*alpha(m,:)')^2;
    cost2 = cost2+temp;

    % constraint c
    Ci = F*alpha(m,:)' >= T(m,:)';
    cons_vec = vertcat(cons_vec,Ci);

    % constraint d
    Di = ( n_min(m) <= p*alpha(m,:)' ) <= n_max(m);
    cons_vec = vertcat(cons_vec,Di);

end

% constraint e

cons_vec = vertcat(cons_vec,M*delta==1);

% objective function

obj = c*cost1 + l*cost2 + 0.5*norm(delta)^2;

% YALMIP

optimize(cons_vec,obj,ops);

d = value(delta);

z = value(alpha);

end
