
%% nominal

Robot1_ID = fopen('./txt_files/Robot1_nom.txt','w');
Robot2_ID = fopen('./txt_files/Robot2_nom.txt','w');
Robot3_ID = fopen('./txt_files/Robot3_nom.txt','w');
Robot4_ID = fopen('./txt_files/Robot4_nom.txt','w');
Robot5_ID = fopen('./txt_files/Robot5_nom.txt','w');
Robot6_ID = fopen('./txt_files/Robot6_nom.txt','w');
Robot7_ID = fopen('./txt_files/Robot7_nom.txt','w');
Robot8_ID = fopen('./txt_files/Robot8_nom.txt','w');
Robot9_ID = fopen('./txt_files/Robot9_nom.txt','w');
Robot10_ID = fopen('./txt_files/Robot10_nom.txt','w');


for i = 1:max(size(tspan_t1))-1

    % task t1
    r5x = x_nom_t1(1,i); r5y = y_nom_t1(1,i);
    r8x = x_nom_t1(2,i); r8y = y_nom_t1(2,i);
    r9x = x_nom_t1(3,i); r9y = y_nom_t1(3,i);

    theta5 = atan2( y_nom_t1(1,i)-y_nom_t1(1,i+1), x_nom_t1(1,i)-x_nom_t1(1,i+1) );
    theta8 = atan2( y_nom_t1(2,i)-y_nom_t1(2,i+1), x_nom_t1(2,i)-x_nom_t1(2,i+1) );
    theta9 = atan2( y_nom_t1(3,i)-y_nom_t1(3,i+1), x_nom_t1(3,i)-x_nom_t1(3,i+1) );

    % task t3
    r2x = x_nom_t2(1,i); r2y = y_nom_t2(1,i);
    r3x = x_nom_t2(2,i); r3y = y_nom_t2(2,i);
    r4x = x_nom_t2(3,i); r4y = y_nom_t2(3,i);
    r6x = x_nom_t2(4,i); r6y = y_nom_t2(4,i);

    theta2 = atan2( y_nom_t2(1,i)-y_nom_t2(1,i+1), x_nom_t2(1,i)-x_nom_t2(1,i+1) );
    theta3 = atan2( y_nom_t2(2,i)-y_nom_t2(2,i+1), x_nom_t2(2,i)-x_nom_t2(2,i+1) );
    theta4 = atan2( y_nom_t2(3,i)-y_nom_t2(3,i+1), x_nom_t2(3,i)-x_nom_t2(3,i+1) );
    theta6 = atan2( y_nom_t2(4,i)-y_nom_t2(4,i+1), x_nom_t2(4,i)-x_nom_t2(4,i+1) );

    % task t3
    r1x = x_nom_t3(1,i); r1y = y_nom_t3(1,i);
    r7x = x_nom_t3(2,i); r7y = y_nom_t3(2,i);
    r10x = x_nom_t3(3,i); r10y = y_nom_t3(3,i);

    theta1 = atan2( y_nom_t3(1,i)-y_nom_t3(1,i+1), x_nom_t3(1,i)-x_nom_t3(1,i+1) );
    theta7 = atan2( y_nom_t3(2,i)-y_nom_t3(2,i+1), x_nom_t3(2,i)-x_nom_t3(2,i+1) );
    theta10 = atan2( y_nom_t3(3,i)-y_nom_t3(3,i+1), x_nom_t3(3,i)-x_nom_t3(3,i+1) );


    fprintf(Robot1_ID,'%f \t %f \t %f \n',r1x,r1y,theta1);
    fprintf(Robot2_ID,'%f \t %f \t %f \n',r2x,r2y,theta2);
    fprintf(Robot3_ID,'%f \t %f \t %f \n',r3x,r3y,theta3);
    fprintf(Robot4_ID,'%f \t %f \t %f \n',r4x,r4y,theta4);
    fprintf(Robot5_ID,'%f \t %f \t %f \n',r5x,r5y,theta5);
    fprintf(Robot6_ID,'%f \t %f \t %f \n',r6x,r6y,theta6);
    fprintf(Robot7_ID,'%f \t %f \t %f \n',r7x,r7y,theta7);
    fprintf(Robot8_ID,'%f \t %f \t %f \n',r8x,r8y,theta8);
    fprintf(Robot9_ID,'%f \t %f \t %f \n',r9x,r9y,theta9);
    fprintf(Robot10_ID,'%f \t %f \t %f \n',r10x,r10y,theta10);

end

fclose(Robot1_ID);
fclose(Robot2_ID);
fclose(Robot3_ID);
fclose(Robot4_ID);
fclose(Robot5_ID);
fclose(Robot6_ID);
fclose(Robot7_ID); 
fclose(Robot8_ID);
fclose(Robot9_ID);
fclose(Robot10_ID); 

%% sat

Robot1_ID = fopen('./txt_files/Robot1_sat.txt','w');
Robot2_ID = fopen('./txt_files/Robot2_sat.txt','w');
Robot3_ID = fopen('./txt_files/Robot3_sat.txt','w');
Robot4_ID = fopen('./txt_files/Robot4_sat.txt','w');
Robot5_ID = fopen('./txt_files/Robot5_sat.txt','w');
Robot6_ID = fopen('./txt_files/Robot6_sat.txt','w');
Robot7_ID = fopen('./txt_files/Robot7_sat.txt','w');
Robot8_ID = fopen('./txt_files/Robot8_sat.txt','w');
Robot9_ID = fopen('./txt_files/Robot9_sat.txt','w');
Robot10_ID = fopen('./txt_files/Robot10_sat.txt','w');


for i = 1:max(size(tspan_t1))-1

    % task t1
    r5x = x_sat_t1(1,i); r5y = y_sat_t1(1,i);
    r8x = x_sat_t1(2,i); r8y = y_sat_t1(2,i);
    r9x = x_sat_t1(3,i); r9y = y_sat_t1(3,i);

    theta5 = atan2( y_sat_t1(1,i)-y_sat_t1(1,i+1), x_sat_t1(1,i)-x_sat_t1(1,i+1) );
    theta8 = atan2( y_sat_t1(2,i)-y_sat_t1(2,i+1), x_sat_t1(2,i)-x_sat_t1(2,i+1) );
    theta9 = atan2( y_sat_t1(3,i)-y_sat_t1(3,i+1), x_sat_t1(3,i)-x_sat_t1(3,i+1) );

    % task t3
    r2x = x_sat_t2(1,i); r2y = y_sat_t2(1,i);
    r3x = x_sat_t2(2,i); r3y = y_sat_t2(2,i);
    r4x = x_sat_t2(3,i); r4y = y_sat_t2(3,i);
    r6x = x_sat_t2(4,i); r6y = y_sat_t2(4,i);

    theta2 = atan2( y_sat_t2(1,i)-y_sat_t2(1,i+1), x_sat_t2(1,i)-x_sat_t2(1,i+1) );
    theta3 = atan2( y_sat_t2(2,i)-y_sat_t2(2,i+1), x_sat_t2(2,i)-x_sat_t2(2,i+1) );
    theta4 = atan2( y_sat_t2(3,i)-y_sat_t2(3,i+1), x_sat_t2(3,i)-x_sat_t2(3,i+1) );
    theta6 = atan2( y_sat_t2(4,i)-y_sat_t2(4,i+1), x_sat_t2(4,i)-x_sat_t2(4,i+1) );

    % task t3
    r1x = x_sat_t3(1,i); r1y = y_sat_t3(1,i);
    r7x = x_sat_t3(2,i); r7y = y_sat_t3(2,i);
    r10x = x_sat_t3(3,i); r10y = y_sat_t3(3,i);

    theta1 = atan2( y_sat_t3(1,i)-y_sat_t3(1,i+1), x_sat_t3(1,i)-x_sat_t3(1,i+1) );
    theta7 = atan2( y_sat_t3(2,i)-y_sat_t3(2,i+1), x_sat_t3(2,i)-x_sat_t3(2,i+1) );
    theta10 = atan2( y_sat_t3(3,i)-y_sat_t3(3,i+1), x_sat_t3(3,i)-x_sat_t3(3,i+1) );


    fprintf(Robot1_ID,'%f \t %f \t %f \n',r1x,r1y,theta1);
    fprintf(Robot2_ID,'%f \t %f \t %f \n',r2x,r2y,theta2);
    fprintf(Robot3_ID,'%f \t %f \t %f \n',r3x,r3y,theta3);
    fprintf(Robot4_ID,'%f \t %f \t %f \n',r4x,r4y,theta4);
    fprintf(Robot5_ID,'%f \t %f \t %f \n',r5x,r5y,theta5);
    fprintf(Robot6_ID,'%f \t %f \t %f \n',r6x,r6y,theta6);
    fprintf(Robot7_ID,'%f \t %f \t %f \n',r7x,r7y,theta7);
    fprintf(Robot8_ID,'%f \t %f \t %f \n',r8x,r8y,theta8);
    fprintf(Robot9_ID,'%f \t %f \t %f \n',r9x,r9y,theta9);
    fprintf(Robot10_ID,'%f \t %f \t %f \n',r10x,r10y,theta10);

end

fclose(Robot1_ID);
fclose(Robot2_ID);
fclose(Robot3_ID);
fclose(Robot4_ID);
fclose(Robot5_ID);
fclose(Robot6_ID);
fclose(Robot7_ID); 
fclose(Robot8_ID);
fclose(Robot9_ID);
fclose(Robot10_ID);

%% dead

Robot1_ID = fopen('./txt_files/Robot1_dead.txt','w');
Robot2_ID = fopen('./txt_files/Robot2_dead.txt','w');
Robot3_ID = fopen('./txt_files/Robot3_dead.txt','w');
Robot4_ID = fopen('./txt_files/Robot4_dead.txt','w');
Robot5_ID = fopen('./txt_files/Robot5_dead.txt','w');
Robot6_ID = fopen('./txt_files/Robot6_dead.txt','w');
Robot7_ID = fopen('./txt_files/Robot7_dead.txt','w');
Robot8_ID = fopen('./txt_files/Robot8_dead.txt','w');
Robot9_ID = fopen('./txt_files/Robot9_dead.txt','w');
Robot10_ID = fopen('./txt_files/Robot10_dead.txt','w');


for i = 1:max(size(tspan_t1))-1

    % task t1
    r5x = x_dead_t1(1,i); r5y = y_dead_t1(1,i);
    r8x = x_dead_t1(2,i); r8y = y_dead_t1(2,i);
    r9x = x_dead_t1(3,i); r9y = y_dead_t1(3,i);

    theta5 = atan2( y_dead_t1(1,i)-y_dead_t1(1,i+1), x_dead_t1(1,i)-x_dead_t1(1,i+1) );
    theta8 = atan2( y_dead_t1(2,i)-y_dead_t1(2,i+1), x_dead_t1(2,i)-x_dead_t1(2,i+1) );
    theta9 = atan2( y_dead_t1(3,i)-y_dead_t1(3,i+1), x_dead_t1(3,i)-x_dead_t1(3,i+1) );

    % task t3
    r2x = x_dead_t2(1,i); r2y = y_dead_t2(1,i);
    r3x = x_dead_t2(2,i); r3y = y_dead_t2(2,i);
    r4x = x_dead_t2(3,i); r4y = y_dead_t2(3,i);
    r6x = x_dead_t2(4,i); r6y = y_dead_t2(4,i);

    theta2 = atan2( y_dead_t2(1,i)-y_dead_t2(1,i+1), x_dead_t2(1,i)-x_dead_t2(1,i+1) );
    theta3 = atan2( y_dead_t2(2,i)-y_dead_t2(2,i+1), x_dead_t2(2,i)-x_dead_t2(2,i+1) );
    theta4 = atan2( y_dead_t2(3,i)-y_dead_t2(3,i+1), x_dead_t2(3,i)-x_dead_t2(3,i+1) );
    theta6 = atan2( y_dead_t2(4,i)-y_dead_t2(4,i+1), x_dead_t2(4,i)-x_dead_t2(4,i+1) );

    % task t3
    r1x = x_dead_t3(1,i); r1y = y_dead_t3(1,i);
    r7x = x_dead_t3(2,i); r7y = y_dead_t3(2,i);
    r10x = x_dead_t3(3,i); r10y = y_dead_t3(3,i);

    theta1 = atan2( y_dead_t3(1,i)-y_dead_t3(1,i+1), x_dead_t3(1,i)-x_dead_t3(1,i+1) );
    theta7 = atan2( y_dead_t3(2,i)-y_dead_t3(2,i+1), x_dead_t3(2,i)-x_dead_t3(2,i+1) );
    theta10 = atan2( y_dead_t3(3,i)-y_dead_t3(3,i+1), x_dead_t3(3,i)-x_dead_t3(3,i+1) );


    fprintf(Robot1_ID,'%f \t %f \t %f \n',r1x,r1y,theta1);
    fprintf(Robot2_ID,'%f \t %f \t %f \n',r2x,r2y,theta2);
    fprintf(Robot3_ID,'%f \t %f \t %f \n',r3x,r3y,theta3);
    fprintf(Robot4_ID,'%f \t %f \t %f \n',r4x,r4y,theta4);
    fprintf(Robot5_ID,'%f \t %f \t %f \n',r5x,r5y,theta5);
    fprintf(Robot6_ID,'%f \t %f \t %f \n',r6x,r6y,theta6);
    fprintf(Robot7_ID,'%f \t %f \t %f \n',r7x,r7y,theta7);
    fprintf(Robot8_ID,'%f \t %f \t %f \n',r8x,r8y,theta8);
    fprintf(Robot9_ID,'%f \t %f \t %f \n',r9x,r9y,theta9);
    fprintf(Robot10_ID,'%f \t %f \t %f \n',r10x,r10y,theta10);

end

fclose(Robot1_ID);
fclose(Robot2_ID);
fclose(Robot3_ID);
fclose(Robot4_ID);
fclose(Robot5_ID);
fclose(Robot6_ID);
fclose(Robot7_ID); 
fclose(Robot8_ID);
fclose(Robot9_ID);
fclose(Robot10_ID);