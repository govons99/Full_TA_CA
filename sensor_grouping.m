function A_new = sensor_grouping(A,ad,n)
%SENSOR_GROUPING Summary of this function goes here
%   Detailed explanation goes here

A_new = A;

[row,col] = size(A);


for h=1:n
    A_temp = A*ad^h;

    for i=1:row
        for j=1:col
            if ( A_new(i,j)==0 && A_temp(i,j)~=0 )
                A_new(i,j) = 1;
            end
        end
    end

    A_new;

end

