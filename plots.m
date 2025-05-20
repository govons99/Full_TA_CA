% plotting stuff

%% control plots

% fisrt task

figure()

tiledlayout(3,2)

%subplot(311)
nexttile
hold on; grid on;
plot(tspan_t1,yc_t1(1:3,:),'lineWidth',2);
title('nominal')
str1 = "$u_{"+num2str(cc{idx1}(1))+"}$";
str2 = "$u_{"+num2str(cc{idx1}(2))+"}$";
str3 = "$u_{"+num2str(cc{idx1}(3))+"}$";
lg = legend(str1,str2,str3, ...
    'numcolumns',3,'interpreter','latex','fontsize',12);
lg.Layout.Tile = 'north';
ylabel("$u_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(tspan_t1,yc_t1(4:end,:),'lineWidth',2);
title('nominal')
ylabel("$u_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(312)
nexttile
hold on; grid on;
plot(tspan_sat_t1,yc_sat_t1(1:3,:),'lineWidth',2);
title('saturation')
ylabel("$u_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(tspan_sat_t1,yc_sat_t1(4:end,:),'lineWidth',2);
title('saturation')
ylabel("$u_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(313)
nexttile
hold on; grid on;
plot(t_dead_t1,u_sat_dead_t1(1:3,:),'lineWidth',2);
title('steady-state opt')
ylabel("$u_x$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(t_dead_t1,u_sat_dead_t1(4:end,:),'lineWidth',2);
title('steady-state opt')
ylabel("$u_y$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

%%saveas(gcf,'./simulation_figures/t1_control.eps','epsc')

% second task

figure()

tiledlayout(3,2)

%subplot(311)
nexttile
hold on; grid on;
plot(tspan_t2,yc_t2(1:4,:),'lineWidth',2);
title('nominal')
str1 = "$u_{"+num2str(cc{idx2}(1))+"}$";
str2 = "$u_{"+num2str(cc{idx2}(2))+"}$";
str3 = "$u_{"+num2str(cc{idx2}(3))+"}$";
str4 = "$u_{"+num2str(cc{idx2}(4))+"}$";
lg = legend(str1,str2,str3,str4 ...
    ,'numcolumns',4,'interpreter','latex','fontsize',12);
lg.Layout.Tile = 'north';
ylabel("$u_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(tspan_t2,yc_t2(5:end,:),'lineWidth',2);
title('nominal')
ylabel("$u_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(312)
nexttile
hold on; grid on;
plot(tspan_sat_t2,yc_sat_t2(1:4,:),'lineWidth',2);
title('saturation')
ylabel("$u_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(tspan_sat_t2,yc_sat_t2(5:end,:),'lineWidth',2);
title('saturation')
ylabel("$u_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(313)
nexttile
hold on; grid on;
plot(t_dead_t2,u_sat_dead_t2(1:4,:),'lineWidth',2);
title('steady-state opt')
ylabel("$u_x$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(t_dead_t2,u_sat_dead_t2(5:end,:),'lineWidth',2);
title('steady-state opt')
ylabel("$u_y$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

%%saveas(gcf,'./simulation_figures/t2_control.eps','epsc')

% third task

figure()

tiledlayout(3,2)

%subplot(311)
nexttile
hold on; grid on;
plot(tspan_t3,yc_t3(1:3,:),'lineWidth',2);
title('nominal')
str1 = "$u_{"+num2str(cc{idx3}(1))+"}$";
str2 = "$u_{"+num2str(cc{idx3}(2))+"}$";
str3 = "$u_{"+num2str(cc{idx3}(3))+"}$";
lg = legend(str1,str2,str3 ...
    ,'numcolumns',3,'interpreter','latex','fontsize',12);
lg.Layout.Tile = 'north';
ylabel("$u_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(tspan_t3,yc_t3(4:end,:),'lineWidth',2);
title('nominal')
ylabel("$u_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(312)
nexttile
hold on; grid on;
plot(tspan_sat_t3,yc_sat_t3(1:3,:),'lineWidth',2);
title('saturation')
ylabel("$u_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(tspan_sat_t3,yc_sat_t3(4:end,:),'lineWidth',2);
title('saturation')
ylabel("$u_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(313)
nexttile
hold on; grid on;
plot(t_dead_t3,u_sat_dead_t3(1:3,:),'lineWidth',2);
title('steady-state opt')
ylabel("$u_x$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

nexttile
hold on; grid on;
plot(t_dead_t3,u_sat_dead_t3(4:end,:),'lineWidth',2);
title('steady-state opt')
ylabel("$u_y$",'interpreter','latex')
xlabel("$Time [s]$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

%%saveas(gcf,'./simulation_figures/t3_control.eps','epsc')

%% error plots

% first task

figure()

tiledlayout(3,2)

% subplot(321)
nexttile
hold on; grid on;
plot(tspan_t1,Lt1*(x_nom_t1-delta_x_t1),'lineWidth',2);
title('nominal')
str1 = "$e_{"+num2str(cc{idx1}(1))+"}$";
str2 = "$e_{"+num2str(cc{idx1}(2))+"}$";
str3 = "$e_{"+num2str(cc{idx1}(3))+"}$";
lg = legend(str1,str2,str3 ...
    ,'numcolumns',min(size(yc_t1)),'interpreter','latex','fontsize',12);
lg.Layout.Tile = 'north';
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(322)
nexttile
hold on; grid on;
plot(tspan_t1,Lt1*(y_nom_t1-delta_y_t1),'lineWidth',2);
title('nominal')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(323)
nexttile
hold on; grid on;
plot(tspan_sat_t1,Lt1*(x_sat_t1-delta_x_t1),'lineWidth',2);
title('saturation')
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(324)
nexttile
hold on; grid on;
plot(tspan_sat_t1,Lt1*(y_sat_t1-delta_y_t1),'lineWidth',2);
title('saturation')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(325)
nexttile
hold on; grid on;
plot(t_dead_t1,Lt1*(x_dead_t1-delta_x_t1),'lineWidth',2);
title('steady-state opt')
xlabel("$Time [s]$",'interpreter','latex')
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(326)
nexttile
hold on; grid on;
plot(t_dead_t1,Lt1*(y_dead_t1-delta_y_t1),'lineWidth',2);
title('steady-state opt')
xlabel("$Time [s]$",'interpreter','latex')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

%%saveas(gcf,'./simulation_figures/t1_error.eps','epsc')

% second task

figure()

tiledlayout(3,2)

% subplot(321)
nexttile
hold on; grid on;
plot(tspan_t2,Lt2*(x_nom_t2-delta_x_t2),'lineWidth',2);
title('nominal')
str1 = "$e_{"+num2str(cc{idx2}(1))+"}$";
str2 = "$e_{"+num2str(cc{idx2}(2))+"}$";
str3 = "$e_{"+num2str(cc{idx2}(3))+"}$";
str4 = "$e_{"+num2str(cc{idx2}(4))+"}$";
lg = legend(str1,str2,str3,str4 ...
    ,'numcolumns',min(size(yc_t2)),'interpreter','latex','fontsize',12);
lg.Layout.Tile = 'north';
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(322)
nexttile
hold on; grid on;
plot(tspan_t2,Lt2*(y_nom_t2-delta_y_t2),'lineWidth',2);
title('nominal')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(323)
nexttile
hold on; grid on;
plot(tspan_sat_t2,Lt2*(x_sat_t2-delta_x_t2),'lineWidth',2);
title('saturation')
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(324)
nexttile
hold on; grid on;
plot(tspan_sat_t2,Lt2*(y_sat_t2-delta_y_t2),'lineWidth',2);
title('saturation')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(325)
nexttile
hold on; grid on;
plot(t_dead_t2,Lt2*(x_dead_t2-delta_x_t2),'lineWidth',2);
title('steady-state opt')
xlabel("$Time [s]$",'interpreter','latex')
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(326)
nexttile
hold on; grid on;
plot(t_dead_t2,Lt2*(y_dead_t2-delta_y_t2),'lineWidth',2);
title('steady-state opt')
xlabel("$Time [s]$",'interpreter','latex')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

%%saveas(gcf,'./simulation_figures/t2_error.eps','epsc')

% third task

figure()

tiledlayout(3,2)

% subplot(321)
nexttile
hold on; grid on;
plot(tspan_t3,Lt3*(x_nom_t3-delta_x_t3),'lineWidth',2);
title('nominal')
str1 = "$e_{"+num2str(cc{idx3}(1))+"}$";
str2 = "$e_{"+num2str(cc{idx3}(2))+"}$";
str3 = "$e_{"+num2str(cc{idx3}(3))+"}$";
lg = legend(str1,str2,str3 ...
    ,'numcolumns',min(size(yc_t3)),'interpreter','latex','fontsize',12);
lg.Layout.Tile = 'north';
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(322)
nexttile
hold on; grid on;
plot(tspan_t3,Lt3*(y_nom_t3-delta_y_t3),'lineWidth',2);
title('nominal')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(323)
nexttile
hold on; grid on;
plot(tspan_sat_t3,Lt3*(x_sat_t3-delta_x_t3),'lineWidth',2);
title('saturation')
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(324)
nexttile
hold on; grid on;
plot(tspan_sat_t3,Lt3*(y_sat_t3-delta_y_t3),'lineWidth',2);
title('saturation')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(325)
nexttile
hold on; grid on;
plot(t_dead_t3,Lt3*(x_dead_t3-delta_x_t3),'lineWidth',2);
title('steady-state opt')
xlabel("$Time [s]$",'interpreter','latex')
ylabel("$e_x$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

% subplot(326)
nexttile
hold on; grid on;
plot(t_dead_t3,Lt3*(y_dead_t3-delta_y_t3),'lineWidth',2);
title('steady-state opt')
xlabel("$Time [s]$",'interpreter','latex')
ylabel("$e_y$",'interpreter','latex')
set(gca,'FontSize',12)
hold off;

%%saveas(gcf,'./simulation_figures/t3_error.eps','epsc')

%% trajectory plots

colors_t1 = [
    [0.00, 0.45, 0.74];   % 1. Blue
    [0.85, 0.33, 0.10];   % 2. Red-Orange
    [0.47, 0.67, 0.19];   % 3. Green
];

colors_t2 = [
    [0.49, 0.18, 0.56];   % 4. Purple
    [0.30, 0.75, 0.93];   % 5. Cyan
    [0.64, 0.08, 0.18];   % 6. Dark Red
    [0.80, 0.40, 0.00];   % 7. Orange
];

colors_t3 = [
    [0.85, 0.65, 0.13];   % 8. Golden Brown
    [0.35, 0.80, 0.45];   % 9. Light Green
    [0.50, 0.10, 0.20];   % 10. Deep Maroon
];

figure()
hold on; grid on;

leg = [];

for i=1:3
    plot(x_nom_t1(i,1),y_nom_t1(i,1),'o','linewidth', 2,'color',colors_t1(i,:));
    plot(x_nom_t1(i,end),y_nom_t1(i,end),'*','linewidth', 2,'color',colors_t1(i,:));
    p1 = plot(x_nom_t1(i,:),y_nom_t1(i,:),'color',colors_t1(i,:));
    leg = vertcat(leg,p1);
end
for i=1:4
    plot(x_nom_t2(i,1),y_nom_t2(i,1),'o','linewidth', 2,'color',colors_t2(i,:));
    plot(x_nom_t2(i,end),y_nom_t2(i,end),'*','linewidth', 2,'color',colors_t2(i,:));
    p1 = plot(x_nom_t2(i,:),y_nom_t2(i,:),'color',colors_t2(i,:));
    leg = vertcat(leg,p1);
end
for i=1:3
    plot(x_nom_t3(i,1),y_nom_t3(i,1),'o','linewidth', 2,'color',colors_t3(i,:));
    plot(x_nom_t3(i,end),y_nom_t3(i,end),'*','linewidth', 2,'color',colors_t3(i,:));
    p1 = plot(x_nom_t3(i,:),y_nom_t3(i,:),'color',colors_t3(i,:));
    leg = vertcat(leg,p1);
end
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
xlim([-12,12])
ylim([-17,10])
%first task
line(x_nom_t1(1:2, end), y_nom_t1(1:2, end),'linewidth', 1.5, 'color', '#a8a8a8');
line(x_nom_t1(2:3, end), y_nom_t1(2:3, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_nom_t1(1, end) x_nom_t1(3, end)], [y_nom_t1(1, end) y_nom_t1(3, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
%second task
line(x_nom_t2(1:2, end), y_nom_t2(1:2, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_nom_t2(1, end) x_nom_t2(3, end)], [y_nom_t2(1, end) y_nom_t2(3, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
line(x_nom_t2(3:4, end), y_nom_t2(3:4, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_nom_t2(2, end) x_nom_t2(4, end)], [y_nom_t2(2, end) y_nom_t2(4, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
axis equal
legend(leg,"$$r_1$$","$$r_2$$","$$r_3$$","$$r_4$$","$$r_5$$","$$r_6$$","$$r_7$$","$$r_8$$","$$r_9$$","$$r_{10}$$", ...
    'location','northoutside','numcolumns',5,'interpreter','latex','fontsize',13);
%saveas(gcf,'./simulation_figures/nom_traj.eps','epsc')

figure()
hold on; grid on;

leg = [];

for i=1:3
    plot(x_sat_t1(i,1),y_sat_t1(i,1),'o','linewidth', 2,'color',colors_t1(i,:));
    plot(x_sat_t1(i,end),y_sat_t1(i,end),'*','linewidth', 2,'color',colors_t1(i,:));
    p1 = plot(x_sat_t1(i,:),y_sat_t1(i,:),'color',colors_t1(i,:));
    leg = vertcat(leg,p1);
end
for i=1:4
    plot(x_sat_t2(i,1),y_sat_t2(i,1),'o','linewidth', 2,'color',colors_t2(i,:));
    plot(x_sat_t2(i,end),y_sat_t2(i,end),'*','linewidth', 2,'color',colors_t2(i,:));
    p1 = plot(x_sat_t2(i,:),y_sat_t2(i,:),'color',colors_t2(i,:));
    leg = vertcat(leg,p1);
end
for i=1:3
    plot(x_sat_t3(i,1),y_sat_t3(i,1),'o','linewidth', 2,'color',colors_t3(i,:));
    plot(x_sat_t3(i,end),y_sat_t3(i,end),'*','linewidth', 2,'color',colors_t3(i,:));
    p1 = plot(x_sat_t3(i,:),y_sat_t3(i,:),'color',colors_t3(i,:));
    leg = vertcat(leg,p1);
end
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
xlim([-12,12])
ylim([-17,10])
%first task
line(x_sat_t1(1:2, end), y_sat_t1(1:2, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line(x_sat_t1(2:3, end), y_sat_t1(2:3, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_sat_t1(1, end) x_sat_t1(3, end)], [y_sat_t1(1, end) y_sat_t1(3, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
%second task
line(x_sat_t2(1:2, end), y_sat_t2(1:2, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_sat_t2(1, end) x_sat_t2(3, end)], [y_sat_t2(1, end) y_sat_t2(3, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
line(x_sat_t2(3:4, end), y_sat_t2(3:4, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_sat_t2(2, end) x_sat_t2(4, end)], [y_sat_t2(2, end) y_sat_t2(4, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
axis equal
legend(leg,"$$r_1$$","$$r_2$$","$$r_3$$","$$r_4$$","$$r_5$$","$$r_6$$","$$r_7$$","$$r_8$$","$$r_9$$","$$r_{10}$$", ...
    'location','northoutside','numcolumns',5,'interpreter','latex','fontsize',13);
%saveas(gcf,'./simulation_figures/sat_traj.eps','epsc')

figure()
hold on; grid on;

leg = [];

for i=1:3
    plot(x_dead_t1(i,1),y_dead_t1(i,1),'o','linewidth', 2,'color',colors_t1(i,:));
    plot(x_dead_t1(i,end),y_dead_t1(i,end),'*','linewidth', 2,'color',colors_t1(i,:));
    p1 = plot(x_dead_t1(i,:),y_dead_t1(i,:),'color',colors_t1(i,:));
    leg = vertcat(leg,p1);
end
for i=1:4
    plot(x_dead_t2(i,1),y_dead_t2(i,1),'o','linewidth', 2,'color',colors_t2(i,:));
    plot(x_dead_t2(i,end),y_dead_t2(i,end),'*','linewidth', 2,'color',colors_t2(i,:));
    p1 = plot(x_dead_t2(i,:),y_dead_t2(i,:),'color',colors_t2(i,:));
    leg = vertcat(leg,p1);
end
for i=1:3
    plot(x_dead_t3(i,1),y_dead_t3(i,1),'o','linewidth', 2,'color',colors_t3(i,:));
    plot(x_dead_t3(i,end),y_dead_t3(i,end),'*','linewidth', 2,'color',colors_t3(i,:));
    p1 = plot(x_dead_t3(i,:),y_dead_t3(i,:),'color',colors_t3(i,:));
    leg = vertcat(leg,p1);
end
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
xlim([-12,12])
ylim([-17,10])
%first task
line(x_dead_t1(1:2, end), y_dead_t1(1:2, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line(x_dead_t1(2:3, end), y_dead_t1(2:3, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_dead_t1(1, end) x_dead_t1(3, end)], [y_dead_t1(1, end) y_dead_t1(3, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
%second task
line(x_dead_t2(1:2, end), y_dead_t2(1:2, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_dead_t2(1, end) x_dead_t2(3, end)], [y_dead_t2(1, end) y_dead_t2(3, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
line(x_dead_t2(3:4, end), y_dead_t2(3:4, end), 'linewidth', 1.5, 'color', '#a8a8a8');
line([x_dead_t2(2, end) x_dead_t2(4, end)], [y_dead_t2(2, end) y_dead_t2(4, end)], 'linewidth', 1.5, 'color', '#a8a8a8');
axis equal
legend(leg,"$$r_1$$","$$r_2$$","$$r_3$$","$$r_4$$","$$r_5$$","$$r_6$$","$$r_7$$","$$r_8$$","$$r_9$$","$$r_{10}$$", ...
    'location','northoutside','numcolumns',5,'interpreter','latex','fontsize',13);
%saveas(gcf,'./simulation_figures/dead_traj.eps','epsc')