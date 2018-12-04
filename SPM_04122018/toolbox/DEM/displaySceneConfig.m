function displaySceneConfig(scene_config_string,lin_length,dot_spread_fac)
% displayAmatrix Displays different scenes using configuration specified by 
% scene_config_string


if sum(strcmp(scene_config_string,'U') ==1) && sum(strcmp(scene_config_string,'R') ==1)
    scene_name = 'UP RIGHT';
elseif sum(strcmp(scene_config_string,'R') ==1) && sum(strcmp(scene_config_string,'D') ==1)
    scene_name = 'RIGHT DOWN';
elseif sum(strcmp(scene_config_string,'D') ==1) && sum(strcmp(scene_config_string,'L') ==1)
    scene_name = 'DOWN LEFT';
elseif sum(strcmp(scene_config_string,'L') ==1) && sum(strcmp(scene_config_string,'U') ==1)
    scene_name = 'LEFT UP';
end


xlimz = [0 1000];
ylimz = [0 1000];

if ~exist('dot_spread_fac','var') || isempty(dot_spread_fac)
    dot_spread_fac = 0.05;
end

dot_spread = [(dot_spread_fac * diff(xlimz)) (dot_spread_fac * diff(ylimz))];


if ~exist('lin_length','var') || isempty(lin_length)
    lin_length = 150;
end

figure('Position',[100 100 1000 700]);
set(gca,'xtick',[])
set(gca,'ytick',[])

hold on;
for loc = 1:9
    
    switch loc
        case 1 %center fixation
            center_coords = [xlimz(2)/2, ylimz(2)/2];
        case 2 % upper left, quadrant 1
            center_coords = [xlimz(2)/4, 3*(ylimz(2)/4)];
        case 3 % lower left, quadrant 2
            center_coords = [xlimz(2)/4, ylimz(2)/4];
        case 4 % upper right, quadrant 3
            center_coords = [3*(xlimz(2)/4), 3*(ylimz(2)/4)];
        case 5 % lower right, quadrant 4
            center_coords = [3*(xlimz(2)/4), ylimz(2)/4];
        case 6 % choose scene A (scene A 'button')
            center_coords = [5*(xlimz(2)/16),ylimz(2)/8];
            text(center_coords(1)-12.5,center_coords(2)-50,'UR','FontSize',18,'FontWeight','bold')
        case 7 % choose scene B (scene B 'button')
            center_coords = [7*(xlimz(2)/16), ylimz(2)/8];
            text(center_coords(1)-12.5,center_coords(2)-50,'RD','FontSize',18,'FontWeight','bold')
        case 8 % choose scene C (scene C 'button')
            center_coords = [9*(xlimz(2)/16), ylimz(2)/8];
            text(center_coords(1)-12.5,center_coords(2)-50,'DL','FontSize',18,'FontWeight','bold')
        case 9 % choose scene D (scene D 'button')
            center_coords = [11*(xlimz(2)/16), ylimz(2)/8];
            text(center_coords(1)-12.5,center_coords(2)-50,'LU','FontSize',18,'FontWeight','bold')
    end
    
    if loc >= 2 && loc <= 5
    
    switch scene_config_string{loc-1}
        case 'n' % 'null'
            scatter(center_coords(1),center_coords(2),50,'ro','filled')
        case 'U' % 'U' == UP == 90 degrees
            
            dots_coords = center_coords + dot_spread.*randn(100,2);
            plot(dots_coords(:,1),dots_coords(:,2),'bo','MarkerSize',5);
            
            angle = 90;
            x(1) = center_coords(1);
            y(1) = center_coords(2);
            x(2) = x(1) + lin_length * cosd(angle);
            y(2) = y(1) + lin_length * sind(angle);
            
            daspect;
            arrow3([x(1) y(1)],[x(2) y(2)],'r-5',3,2,[],0.75);
            
        case 'R' % 'R' == RIGHT == 0 degrees
            
            dots_coords = center_coords + dot_spread.*randn(100,2);
            plot(dots_coords(:,1),dots_coords(:,2),'bo','MarkerSize',5);
            
            angle = 0;
            x(1) = center_coords(1);
            y(1) = center_coords(2);
            x(2) = x(1) + lin_length * cosd(angle);
            y(2) = y(1) + lin_length * sind(angle);
            
            daspect;
            arrow3([x(1) y(1)],[x(2) y(2)],'r-5',3,2,[],0.75);
            
           
        case 'D' % 'R' == DOWN == 270 degrees
            
            dots_coords = center_coords + dot_spread.*randn(100,2);
            plot(dots_coords(:,1),dots_coords(:,2),'bo','MarkerSize',5);
            
            angle = 270;
            x(1) = center_coords(1);
            y(1) = center_coords(2);
            x(2) = x(1) + lin_length * cosd(angle);
            y(2) = y(1) + lin_length * sind(angle);
            
            daspect;
            arrow3([x(1) y(1)],[x(2) y(2)],'r-5',3,2,[],0.75);
            
        case 'L' % 'L' == LEFT == 180 degrees
            
            dots_coords = center_coords + dot_spread.*randn(100,2);
            plot(dots_coords(:,1),dots_coords(:,2),'bo','MarkerSize',5);
            
            angle = 180;
            x(1) = center_coords(1);
            y(1) = center_coords(2);
            x(2) = x(1) + lin_length * cosd(angle);
            y(2) = y(1) + lin_length * sind(angle);
            
            daspect;
            arrow3([x(1) y(1)],[x(2) y(2)],'r-5',3,2,[],0.75);
            
        case 6 % correct
            scatter(center_coords(1),center_coords(2),175,'go','filled')
        case 7 % incorrect
            scatter(center_coords(1),center_coords(2),175,'ro')
    end
    xlim(xlimz)
    ylim(ylimz);
end

horz_line = [ xlimz(1) xlimz(end) ; ylimz(end)/2 ylimz(end)/2];
vert_line = [ xlimz(end)/2 xlimz(end)/2; ylimz(1) ylimz(end) ];

plot(horz_line(1,:),horz_line(2,:),'k--','LineWidth',1.5); 
plot(vert_line(1,:),vert_line(2,:),'k--','LineWidth',1.5);

title(sprintf('Category: %s',scene_name))


end