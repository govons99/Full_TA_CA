function S = specialization(T,F,nt)

S = zeros(nt,nt);

for i=1:nt
    
    if ( T(i,:)*F > 0 )
        S(i,i) = 1;
    end

end

end