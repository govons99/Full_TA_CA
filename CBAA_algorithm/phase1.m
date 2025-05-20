function [x,y,Ji] = phase1(c,x_old,y_old,J)

    Ji = J;
    
    % number of tasks
    nt = max(size(c));
    
    % list of valid tasks
    hi = zeros(nt,1);
    
    % agent's task list
    x = x_old;

    % winning bids list
    y = y_old;
    
    is_assigned = sum(x);

    % check if i-th agent is unassigned
    if is_assigned==0
        % create the valid task list
        for j=1:nt
            if (c(j)>y(j))
                hi(j) = 1;
            end
        end

        if any(hi~=0)
            % select the task based on the winning bids
            [~,Ji] = max(hi.*c);
            x(Ji) = 1;
            y(Ji) = c(Ji);
        end
    else
        disp('already assigned')
    end


end
