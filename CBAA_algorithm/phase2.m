function [xi,yi] = phase2(x,yi,Y,Ji,idx)
    
    xi = x;

    [~,nt] = size(Y);

    for j =1:nt
        yi(j) = max(Y(:,j));
    end

    [bar_z,z_iJ] = max(Y(:,Ji));
    
    Y
    Ji

    if yi(Ji) == bar_z
        xi(Ji) = 1;
    else
        xi(Ji) = 0;
    end


    % if Ji==1
    %     [~,z_iJ] = max(Y(:,Ji).*r1);
    % elseif Ji==2
    %     [~,z_iJ] = max(Y(:,Ji).*r2);
    % else
    %     [~,z_iJ] = max(Y(:,Ji).*r3);
    % end
    
    % if z_iJ~=idx
    %     xi(Ji) = 0;
    % end
    

end