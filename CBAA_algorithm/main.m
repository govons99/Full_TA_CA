clc
clear
close all

% agent 1
x1_ini = [0, 0, 0];
y1_ini = [0, 0, 0];
c1 = [8; 0; 1];
J1_ini = 0;

% agent 2
x2_ini = [0, 0, 0];
y2_ini = [0, 0, 0];
c2 = [4; 8; 0];
J2_ini = 0;

% agent 3
x3_ini = [0, 0, 0];
y3_ini = [0, 0, 0];
%c3 = [1; 1; 10];
c3 = c1;
J3_ini = 0;

for j = 1:2

    % phase1
    
    [x1_1,y1_1,J1] = phase1(c1,x1_ini,y1_ini,J1_ini)
    [x2_1,y2_1,J2] = phase1(c2,x2_ini,y2_ini,J2_ini)
    [x3_1,y3_1,J3] = phase1(c3,x3_ini,y3_ini,J3_ini)
    
    % phase2
    
    [x1_2,y1_2] = phase2(x1_1,y1_1,[y1_1;y2_1;y3_1],J1,1)
    [x2_2,y2_2] = phase2(x2_1,y2_1,[y1_1;y2_1;y3_1],J2,2)
    [x3_2,y3_2] = phase2(x3_1,y3_1,[y1_1;y2_1;y3_1],J3,3)

    x1_ini = x1_2; y1_ini = y1_2;
    x2_ini = x2_2; y2_ini = y2_2;
    x3_ini = x3_2; y3_ini = y3_2;

    J1_ini = J1;
    J2_ini = J2;
    J3_ini = J3;

end
