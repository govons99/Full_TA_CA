function F = Rob2Cap_Matrix(H,A,nc,nr)

F = zeros(nc,nr);

for i=1:max(size(H))
    fi = H{i}*A;
    F(i,:) = max_kron(1,fi);
end

end