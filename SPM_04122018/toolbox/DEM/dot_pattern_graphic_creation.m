%% create dot pattern graphics for visualization of dot pattern trials

center = [0 0];
radius = 0.5;
ang = 0:0.1:(2*pi);

boundary_x = center(1) + radius*cos(ang);
boundary_y = center(2) + radius*sin(ang);

num_dots = 500;
low_edge = -0.5;
high_edge = 0.5;

% CREATE THE 'UP' PATTERN
dot_coords = low_edge + (high_edge - low_edge).* rand(500,2);
dot_coords = dot_coords(inpolygon(dot_coords(:,1),dot_coords(:,2),boundary_x,boundary_y),:);
figure('Position',[100 800 350 300])
scatter(dot_coords(:,1),dot_coords(:,2),200,'k.')
hold on; arrow3([0 -0.3],[0 0.3],'b-7',7,5,[],1);
set(gca,'Visible','off')
currIm = getframe(gca);
UP = currIm.cdata;
close gcf;

% CREATE THE 'RIGHT' PATTERN
dot_coords = low_edge + (high_edge - low_edge).* rand(500,2);
dot_coords = dot_coords(inpolygon(dot_coords(:,1),dot_coords(:,2),boundary_x,boundary_y),:);
figure('Position',[100 800 350 300])
scatter(dot_coords(:,1),dot_coords(:,2),200,'k.')
hold on; arrow3([-0.3 0],[0.3 0],'b-7',7,5,[],1);
set(gca,'Visible','off')
currIm = getframe(gca);
RIGHT = currIm.cdata;
close gcf;

% CREATE THE 'DOWN' PATTERN
dot_coords = low_edge + (high_edge - low_edge).* rand(500,2);
dot_coords = dot_coords(inpolygon(dot_coords(:,1),dot_coords(:,2),boundary_x,boundary_y),:);
figure('Position',[100 800 350 300])
scatter(dot_coords(:,1),dot_coords(:,2),200,'k.')
hold on; arrow3([0 0.3],[0 -0.3],'b-7',7,5,[],1);
set(gca,'Visible','off')
currIm = getframe(gca);
DOWN = currIm.cdata;
close gcf;

% CREATE THE 'LEFT' PATTERN
dot_coords = low_edge + (high_edge - low_edge).* rand(500,2);
dot_coords = dot_coords(inpolygon(dot_coords(:,1),dot_coords(:,2),boundary_x,boundary_y),:);
figure('Position',[100 800 350 300])
scatter(dot_coords(:,1),dot_coords(:,2),200,'k.')
hold on; arrow3([0.3 0],[-0.3 0],'b-7',7,5,[],1);
set(gca,'Visible','off')
currIm = getframe(gca);
LEFT = currIm.cdata;
close gcf;

save('MDP_search_graphics_DOT.mat','UP','RIGHT','DOWN','LEFT')