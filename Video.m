close all

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

%% Nominal

video_nom = VideoWriter('nominal_traj.avi');
video_nom.Quality = 95;
open(video_nom);

figure()
hold on; grid on;

xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
xlim([-12,12])
ylim([-17,10])
axis equal

leg = [];
star = [];

% plotting initial pos
for i=1:3
    plot(x_nom_t1(i,1),y_nom_t1(i,1),'o','linewidth', 2,'color',colors_t1(i,:));
    plot(x_nom_t2(i,1),y_nom_t2(i,1),'o','linewidth', 2,'color',colors_t2(i,:));
    plot(x_nom_t3(i,1),y_nom_t3(i,1),'o','linewidth', 2,'color',colors_t3(i,:));
end
plot(x_nom_t2(4,1),y_nom_t2(4,1),'o','linewidth', 2,'color',colors_t2(4,:));

% plotting the evolution

dim = max(size(tspan_t1));

j = 1;

for i=1:3

    % current position
    s1 = plot(x_nom_t1(i,j),y_nom_t1(i,j),'*','linewidth', 2,'color',colors_t1(i,:));
    s2 = plot(x_nom_t2(i,j),y_nom_t2(i,j),'*','linewidth', 2,'color',colors_t2(i,:));
    s3 = plot(x_nom_t3(i,j),y_nom_t3(i,j),'*','linewidth', 2,'color',colors_t3(i,:));
    star = vertcat(star,[s1;s2;s3]);

    % trajectory
    p1 = plot(x_nom_t1(i,1:j),y_nom_t1(i,1:j),'color',colors_t1(i,:));
    p2 = plot(x_nom_t2(i,1:j),y_nom_t2(i,1:j),'color',colors_t2(i,:));
    p3 = plot(x_nom_t3(i,1:j),y_nom_t3(i,1:j),'color',colors_t3(i,:));
    leg = vertcat(leg,[p1;p2;p3]);

end
s4 = plot(x_nom_t2(4,j),y_nom_t2(4,j),'*','linewidth', 2,'color',colors_t2(4,:));
star = vertcat(star,s4);
p4 = plot(x_nom_t2(4,1:j),y_nom_t2(4,1:j),'color',colors_t2(4,:));
leg = vertcat(leg,p4);

%first task
l1 = line(x_nom_t1(1:2, j), y_nom_t1(1:2, j),'linewidth', 1.5, 'color', '#a8a8a8');
l2 = line(x_nom_t1(2:3, j), y_nom_t1(2:3, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l3 = line([x_nom_t1(1, j) x_nom_t1(3, j)], [y_nom_t1(1, j) y_nom_t1(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
%second task
l4 = line(x_nom_t2(1:2, j), y_nom_t2(1:2, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l5 = line([x_nom_t2(1, j) x_nom_t2(3, j)], [y_nom_t2(1, j) y_nom_t2(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
l6 = line(x_nom_t2(3:4, j), y_nom_t2(3:4, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l7 = line([x_nom_t2(2, j) x_nom_t2(4, j)], [y_nom_t2(2, j) y_nom_t2(4, j)], 'linewidth', 1.5, 'color', '#a8a8a8');

for j=1:20:dim

    delete(l1)
    delete(l2)
    delete(l3)
    delete(l4)
    delete(l5)
    delete(l6)
    delete(l7)
    delete(star)

    
    for i=1:3

        % current position
        s1 = plot(x_nom_t1(i,j),y_nom_t1(i,j),'*','linewidth', 2,'color',colors_t1(i,:));
        s2 = plot(x_nom_t2(i,j),y_nom_t2(i,j),'*','linewidth', 2,'color',colors_t2(i,:));
        s3 = plot(x_nom_t3(i,j),y_nom_t3(i,j),'*','linewidth', 2,'color',colors_t3(i,:));
        star = vertcat(star,[s1;s2;s3]);

        % trajectory
        plot(x_nom_t1(i,1:j),y_nom_t1(i,1:j),'color',colors_t1(i,:));
        plot(x_nom_t2(i,1:j),y_nom_t2(i,1:j),'color',colors_t2(i,:));
        plot(x_nom_t3(i,1:j),y_nom_t3(i,1:j),'color',colors_t3(i,:));

    end
    s4 = plot(x_nom_t2(4,j),y_nom_t2(4,j),'*','linewidth', 2,'color',colors_t2(4,:));
    star = vertcat(star,s4);
    plot(x_nom_t2(4,1:j),y_nom_t2(4,1:j),'color',colors_t2(4,:));

    %first task
    l1 = line(x_nom_t1(1:2, j), y_nom_t1(1:2, j),'linewidth', 1.5, 'color', '#a8a8a8');
    l2 = line(x_nom_t1(2:3, j), y_nom_t1(2:3, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l3 = line([x_nom_t1(1, j) x_nom_t1(3, j)], [y_nom_t1(1, j) y_nom_t1(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
    %second task
    l4 = line(x_nom_t2(1:2, j), y_nom_t2(1:2, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l5 = line([x_nom_t2(1, j) x_nom_t2(3, j)], [y_nom_t2(1, j) y_nom_t2(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
    l6 = line(x_nom_t2(3:4, j), y_nom_t2(3:4, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l7 = line([x_nom_t2(2, j) x_nom_t2(4, j)], [y_nom_t2(2, j) y_nom_t2(4, j)], 'linewidth', 1.5, 'color', '#a8a8a8');

    legend(leg,"$$r_1$$","$$r_2$$","$$r_3$$","$$r_4$$","$$r_5$$","$$r_6$$","$$r_7$$","$$r_8$$","$$r_9$$","$$r_{10}$$", ...
    'location','northoutside','numcolumns',5,'interpreter','latex','fontsize',13);

    frame = getframe(gcf);
    writeVideo(video_nom,frame);

    drawnow

end

hold off

frame = getframe(gcf);
writeVideo(video_nom,frame);

close(video_nom);

%% Saturating

video_sat = VideoWriter('saturation_traj.avi');
video_sat.Quality = 95;
open(video_sat);

figure()
hold on; grid on;

xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
xlim([-12,12])
ylim([-17,10])
axis equal

leg = [];
star = [];

% plotting initial pos
for i=1:3
    plot(x_sat_t1(i,1),y_sat_t1(i,1),'o','linewidth', 2,'color',colors_t1(i,:));
    plot(x_sat_t2(i,1),y_sat_t2(i,1),'o','linewidth', 2,'color',colors_t2(i,:));
    plot(x_sat_t3(i,1),y_sat_t3(i,1),'o','linewidth', 2,'color',colors_t3(i,:));
end
plot(x_sat_t2(4,1),y_sat_t2(4,1),'o','linewidth', 2,'color',colors_t2(4,:));

% plotting the evolution

dim = max(size(tspan_t1));

j = 1;

for i=1:3

    % current position
    s1 = plot(x_sat_t1(i,j),y_sat_t1(i,j),'*','linewidth', 2,'color',colors_t1(i,:));
    s2 = plot(x_sat_t2(i,j),y_sat_t2(i,j),'*','linewidth', 2,'color',colors_t2(i,:));
    s3 = plot(x_sat_t3(i,j),y_sat_t3(i,j),'*','linewidth', 2,'color',colors_t3(i,:));
    star = vertcat(star,[s1;s2;s3]);

    % trajectory
    p1 = plot(x_sat_t1(i,1:j),y_sat_t1(i,1:j),'color',colors_t1(i,:));
    p2 = plot(x_sat_t2(i,1:j),y_sat_t2(i,1:j),'color',colors_t2(i,:));
    p3 = plot(x_sat_t3(i,1:j),y_sat_t3(i,1:j),'color',colors_t3(i,:));
    leg = vertcat(leg,[p1;p2;p3]);

end
s4 = plot(x_sat_t2(4,j),y_sat_t2(4,j),'*','linewidth', 2,'color',colors_t2(4,:));
star = vertcat(star,s4);
p4 = plot(x_sat_t2(4,1:j),y_sat_t2(4,1:j),'color',colors_t2(4,:));
leg = vertcat(leg,p4);

%first task
l1 = line(x_sat_t1(1:2, j), y_sat_t1(1:2, j),'linewidth', 1.5, 'color', '#a8a8a8');
l2 = line(x_sat_t1(2:3, j), y_sat_t1(2:3, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l3 = line([x_sat_t1(1, j) x_sat_t1(3, j)], [y_sat_t1(1, j) y_sat_t1(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
%second task
l4 = line(x_sat_t2(1:2, j), y_sat_t2(1:2, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l5 = line([x_sat_t2(1, j) x_sat_t2(3, j)], [y_sat_t2(1, j) y_sat_t2(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
l6 = line(x_sat_t2(3:4, j), y_sat_t2(3:4, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l7 = line([x_sat_t2(2, j) x_sat_t2(4, j)], [y_sat_t2(2, j) y_sat_t2(4, j)], 'linewidth', 1.5, 'color', '#a8a8a8');

for j=1:20:dim

    delete(l1)
    delete(l2)
    delete(l3)
    delete(l4)
    delete(l5)
    delete(l6)
    delete(l7)
    delete(star)

    
    for i=1:3

        % current position
        s1 = plot(x_sat_t1(i,j),y_sat_t1(i,j),'*','linewidth', 2,'color',colors_t1(i,:));
        s2 = plot(x_sat_t2(i,j),y_sat_t2(i,j),'*','linewidth', 2,'color',colors_t2(i,:));
        s3 = plot(x_sat_t3(i,j),y_sat_t3(i,j),'*','linewidth', 2,'color',colors_t3(i,:));
        star = vertcat(star,[s1;s2;s3]);

        % trajectory
        plot(x_sat_t1(i,1:j),y_sat_t1(i,1:j),'color',colors_t1(i,:));
        plot(x_sat_t2(i,1:j),y_sat_t2(i,1:j),'color',colors_t2(i,:));
        plot(x_sat_t3(i,1:j),y_sat_t3(i,1:j),'color',colors_t3(i,:));

    end
    s4 = plot(x_sat_t2(4,j),y_sat_t2(4,j),'*','linewidth', 2,'color',colors_t2(4,:));
    star = vertcat(star,s4);
    plot(x_sat_t2(4,1:j),y_sat_t2(4,1:j),'color',colors_t2(4,:));

    %first task
    l1 = line(x_sat_t1(1:2, j), y_sat_t1(1:2, j),'linewidth', 1.5, 'color', '#a8a8a8');
    l2 = line(x_sat_t1(2:3, j), y_sat_t1(2:3, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l3 = line([x_sat_t1(1, j) x_sat_t1(3, j)], [y_sat_t1(1, j) y_sat_t1(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
    %second task
    l4 = line(x_sat_t2(1:2, j), y_sat_t2(1:2, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l5 = line([x_sat_t2(1, j) x_sat_t2(3, j)], [y_sat_t2(1, j) y_sat_t2(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
    l6 = line(x_sat_t2(3:4, j), y_sat_t2(3:4, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l7 = line([x_sat_t2(2, j) x_sat_t2(4, j)], [y_sat_t2(2, j) y_sat_t2(4, j)], 'linewidth', 1.5, 'color', '#a8a8a8');

    legend(leg,"$$r_1$$","$$r_2$$","$$r_3$$","$$r_4$$","$$r_5$$","$$r_6$$","$$r_7$$","$$r_8$$","$$r_9$$","$$r_{10}$$", ...
    'location','northoutside','numcolumns',5,'interpreter','latex','fontsize',13);

    frame = getframe(gcf);
    writeVideo(video_sat,frame);

    drawnow

end

hold off

frame = getframe(gcf);
writeVideo(video_sat,frame);

close(video_sat);


%% Deadzone

video_dead = VideoWriter('deadzone_traj.avi');
video_dead.Quality = 95;
open(video_dead);

figure()
hold on; grid on;

xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
xlim([-12,12])
ylim([-17,10])
axis equal

leg = [];
star = [];

% plotting initial pos
for i=1:3
    plot(x_dead_t1(i,1),y_dead_t1(i,1),'o','linewidth', 2,'color',colors_t1(i,:));
    plot(x_dead_t2(i,1),y_dead_t2(i,1),'o','linewidth', 2,'color',colors_t2(i,:));
    plot(x_dead_t3(i,1),y_dead_t3(i,1),'o','linewidth', 2,'color',colors_t3(i,:));
end
plot(x_dead_t2(4,1),y_dead_t2(4,1),'o','linewidth', 2,'color',colors_t2(4,:));

% plotting the evolution

dim = max(size(tspan_t1));

j = 1;

for i=1:3

    % current position
    s1 = plot(x_dead_t1(i,j),y_dead_t1(i,j),'*','linewidth', 2,'color',colors_t1(i,:));
    s2 = plot(x_dead_t2(i,j),y_dead_t2(i,j),'*','linewidth', 2,'color',colors_t2(i,:));
    s3 = plot(x_dead_t3(i,j),y_dead_t3(i,j),'*','linewidth', 2,'color',colors_t3(i,:));
    star = vertcat(star,[s1;s2;s3]);

    % trajectory
    p1 = plot(x_dead_t1(i,1:j),y_dead_t1(i,1:j),'color',colors_t1(i,:));
    p2 = plot(x_dead_t2(i,1:j),y_dead_t2(i,1:j),'color',colors_t2(i,:));
    p3 = plot(x_dead_t3(i,1:j),y_dead_t3(i,1:j),'color',colors_t3(i,:));
    leg = vertcat(leg,[p1;p2;p3]);

end
s4 = plot(x_dead_t2(4,j),y_dead_t2(4,j),'*','linewidth', 2,'color',colors_t2(4,:));
star = vertcat(star,s4);
p4 = plot(x_dead_t2(4,1:j),y_dead_t2(4,1:j),'color',colors_t2(4,:));
leg = vertcat(leg,p4);

%first task
l1 = line(x_dead_t1(1:2, j), y_dead_t1(1:2, j),'linewidth', 1.5, 'color', '#a8a8a8');
l2 = line(x_dead_t1(2:3, j), y_dead_t1(2:3, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l3 = line([x_dead_t1(1, j) x_dead_t1(3, j)], [y_dead_t1(1, j) y_dead_t1(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
%second task
l4 = line(x_dead_t2(1:2, j), y_dead_t2(1:2, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l5 = line([x_dead_t2(1, j) x_dead_t2(3, j)], [y_dead_t2(1, j) y_dead_t2(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
l6 = line(x_dead_t2(3:4, j), y_dead_t2(3:4, j), 'linewidth', 1.5, 'color', '#a8a8a8');
l7 = line([x_dead_t2(2, j) x_dead_t2(4, j)], [y_dead_t2(2, j) y_dead_t2(4, j)], 'linewidth', 1.5, 'color', '#a8a8a8');

for j=1:20:dim

    delete(l1)
    delete(l2)
    delete(l3)
    delete(l4)
    delete(l5)
    delete(l6)
    delete(l7)
    delete(star)

    
    for i=1:3

        % current position
        s1 = plot(x_dead_t1(i,j),y_dead_t1(i,j),'*','linewidth', 2,'color',colors_t1(i,:));
        s2 = plot(x_dead_t2(i,j),y_dead_t2(i,j),'*','linewidth', 2,'color',colors_t2(i,:));
        s3 = plot(x_dead_t3(i,j),y_dead_t3(i,j),'*','linewidth', 2,'color',colors_t3(i,:));
        star = vertcat(star,[s1;s2;s3]);

        % trajectory
        plot(x_dead_t1(i,1:j),y_dead_t1(i,1:j),'color',colors_t1(i,:));
        plot(x_dead_t2(i,1:j),y_dead_t2(i,1:j),'color',colors_t2(i,:));
        plot(x_dead_t3(i,1:j),y_dead_t3(i,1:j),'color',colors_t3(i,:));

    end
    s4 = plot(x_dead_t2(4,j),y_dead_t2(4,j),'*','linewidth', 2,'color',colors_t2(4,:));
    star = vertcat(star,s4);
    plot(x_dead_t2(4,1:j),y_dead_t2(4,1:j),'color',colors_t2(4,:));

    %first task
    l1 = line(x_dead_t1(1:2, j), y_dead_t1(1:2, j),'linewidth', 1.5, 'color', '#a8a8a8');
    l2 = line(x_dead_t1(2:3, j), y_dead_t1(2:3, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l3 = line([x_dead_t1(1, j) x_dead_t1(3, j)], [y_dead_t1(1, j) y_dead_t1(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
    %second task
    l4 = line(x_dead_t2(1:2, j), y_dead_t2(1:2, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l5 = line([x_dead_t2(1, j) x_dead_t2(3, j)], [y_dead_t2(1, j) y_dead_t2(3, j)], 'linewidth', 1.5, 'color', '#a8a8a8');
    l6 = line(x_dead_t2(3:4, j), y_dead_t2(3:4, j), 'linewidth', 1.5, 'color', '#a8a8a8');
    l7 = line([x_dead_t2(2, j) x_dead_t2(4, j)], [y_dead_t2(2, j) y_dead_t2(4, j)], 'linewidth', 1.5, 'color', '#a8a8a8');

    legend(leg,"$$r_1$$","$$r_2$$","$$r_3$$","$$r_4$$","$$r_5$$","$$r_6$$","$$r_7$$","$$r_8$$","$$r_9$$","$$r_{10}$$", ...
    'location','northoutside','numcolumns',5,'interpreter','latex','fontsize',13);

    frame = getframe(gcf);
    writeVideo(video_dead,frame);

    drawnow

end

hold off

frame = getframe(gcf);
writeVideo(video_dead,frame);

close(video_dead);