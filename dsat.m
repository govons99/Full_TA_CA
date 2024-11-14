function ds = dsat(u,m,r)
% ds=derivative of s
% m = saturation level, r = boundary layer 
if abs(u)>m
   ds=0;
else if abs(u)<m-r
   ds=1;
else
a=sign(u)*2/r^3;
b=-3*(2*m-r)/r^3;
c=sign(u)*6*(m^2-m*r)/r^3;
d=-(2*m^3-3*m^2*r)/r^3;
ds=a.*u.^3+b.*u.^2+c.*u+d;
end

end