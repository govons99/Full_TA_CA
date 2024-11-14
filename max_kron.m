function y = max_kron(n,x)

[row,col] = size(x);

temp = zeros(row,col);

y = zeros(1,col);

for i = 1:row    
    for j = 1:col       
        if ( x(i,j)==n )            
            temp(i,j) = 1;            
        end        
    end
end

for j = 1:col  
    max = 0;
    for i =1:row
        if ( temp(i,j) >= max )
            max = temp(i,j);
        end
    end
    y(j) = max;
end

end