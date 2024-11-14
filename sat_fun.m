function u_sat = sat_fun(u,u_max)
    
    u_sat = u;
    for i=1:max(size(u))
        if (abs(u(i))>u_max(i))
            u_sat(i) = sign(u(i))*u_max(i);
        end
    end

end