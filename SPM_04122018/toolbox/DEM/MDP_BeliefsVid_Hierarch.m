function [ M,dot_idx ] = MDP_BeliefsVid_Hierarch(MDP,trial_idx)
% MDP_BELIEFSVID_HIERARCH 
% Creates video of the trajectory of an MDP active inference agent's
% visual search behavior during a single trial, as well as the
% corresponding evolution of posterior beliefs (at multiple hierarchical
% levels)
% created February 6 2019

all_configs = generate_scenes();
dot_idx = [];

nDots = 20;
apertureSize = [3.5 3.5];
lifetime = 12;

mdp = MDP(1,trial_idx);

outcomes = {'null','UP','RIGHT','DOWN','LEFT'};

scene_config    = all_configs(mdp.s(1,1),:);
scene_config(strcmp(scene_config,'n')) = outcomes(1);
scene_config(strcmp(scene_config,'U')) = outcomes(2);
scene_config(strcmp(scene_config,'R')) = outcomes(3);
scene_config(strcmp(scene_config,'D')) = outcomes(4);
scene_config(strcmp(scene_config,'L')) = outcomes(5);

% load images
%----------------------------------------------------------------------
load MDP_search_graphics_DOT_newest

UP = double(UP);
RIGHT = double(RIGHT);
DOWN = double(DOWN);
LEFT = double(LEFT);
null = zeros(size(UP)) + 1;

% choices: this text will be used to label scene categories at the bottom
%----------------------------------------------------------------------
choice = {'UR','RD','DL','LU'};

T_Higher = find(mdp.o(3,:) > 5,1);
if isempty(T_Higher)
    T_Higher = size(mdp.o,2);
end

%% Plot two subplots within one larger figure

figure('Position',[200 700 1100 750]);

% subplot 1: beliefs about hidden scenes
subplot(1,2,1);

% subplot 2: visual outcomes (the sampled, noisy dot patterns)
subplot(1,2,2);


subplot(1,2,1);
axis([-8 8 -12 8])

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

quadrant_centers = [-5 5;   % centers of the quadrants where the current posterior beliefs will be plotted
                    -5 -5;
                    5 5;
                    5 -5];
text_coords = [-5  -10;
               -1.75  -10;
                1.75  -10;
                5  -10];

%arrays that store vertices for RDM pattern-quadrants in second subplot
quadrantsX = zeros(4,size(quadrant_centers,1)); 
quadrantsY = zeros(4,size(quadrant_centers,1));
    
for patt_i = 1:4
    
    quadrantsX(:,patt_i) = [quadrant_centers(patt_i,1) - 2.5; quadrant_centers(patt_i,1) - 2.5; ...
                            quadrant_centers(patt_i,1) + 2.5; quadrant_centers(patt_i,1) + 2.5];                           
    quadrantsY(:,patt_i) = [quadrant_centers(patt_i,2) - 2.5; quadrant_centers(patt_i,2) + 2.5;...
                            quadrant_centers(patt_i,2) + 2.5; quadrant_centers(patt_i,2) - 2.5];
                        
end

subplot(122);  
patch(quadrantsX,quadrantsY,'black')
axis image xy
axis([-8 8 -12 8])


% x/y coordinates of saccade-able points (starting with fixation center and
% going to scene-categorization choices) 

all_fixatable_points = [ [0,0] ; quadrant_centers ; [text_coords(:,1),text_coords(:,2) + 1] ];

sacc_seq = all_fixatable_points(mdp.o(3,:),:);

for j = 1:2
    T(:,j) = interp(sacc_seq(:,j),8,2);
    T(:,j) = T(:,j) + spm_conv(0.8*randn(size(T(:,j))),2)/16; % doesn't do the eye-movement biomechanics convolution for categorization choices
end

choose_idx = find(mdp.o(3,:) >= 6);

if ~isempty(choose_idx)
    T(choose_idx(1)*8,:) = sacc_seq(choose_idx(1),:);   
    for choice_i = 1:length(choose_idx)-1
        idx = (choose_idx(choice_i)*8 + 1) : (choose_idx(choice_i)*8 + 8);
%         T( (choose_idx(choice_i)*8 + 1):(choose_idx(choice_i+1)*8),:) = repmat(sacc_seq(choose_idx(choice_i),:),8,1);
        T(idx,:) = repmat(sacc_seq(choose_idx(choice_i),:),8,1);
    end
end

quadrant_arrival_times = [1,8:8:length(T)];

%%

frame_num = 1;

for t_high_i = 1:T_Higher
    
    
    subplot(121)
    title(sprintf('Timestep # %d',t_high_i))
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    subplot(122)
    title(sprintf('Timestep # %d',t_high_i))
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    
    for j = 1:4
        S{j} = zeros(size(UP));
    end
    
    for config_i = 1:size(all_configs,1)
        
        a = all_configs(config_i,:);
        a(strcmp(a,'n')) = outcomes(1);
        a(strcmp(a,'U')) = outcomes(2);
        a(strcmp(a,'R')) = outcomes(3);
        a(strcmp(a,'D')) = outcomes(4);
        a(strcmp(a,'L')) = outcomes(5);
        
        % posterior beliefs about scene
        %--------------------------------------------------
        p     = mdp.xn{1}(end,config_i,t_high_i,t_high_i);
        
        for j = 1:4
            S{j} = S{j} + p.*eval(a{j});
        end
        
    end
    
    %% plot current posterior beliefs about hidden states at higher level (scene beliefs) 
    
    subplot(121)
    
    for j = 1:4
        imagesc([quadrant_centers(j,1) - 2.5, quadrant_centers(j,1) + 2.5],...
              [quadrant_centers(j,2) - 2.5, quadrant_centers(j,2) + 2.5], flipud(S{j}/max(S{j}(:)))); % have to flip because of call 'axis image xy' (see a line 137 below)
        hold on;
    end
   
    axis image xy
    axis([-8 8 -12 8])
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    for j = 1:4
        if ismember(mdp.s(1,1),(12 * (j-1) +1): (12*j))
            text(text_coords(j,1),text_coords(j,2),choice{j},'FontSize',14,'FontWeight','Bold','color','blue')
        else
            text(text_coords(j,1),text_coords(j,2),choice{j},'FontSize',14)
        end
    end
    
    
    subplot(122);
    
    for j = 1:4
        if ismember(mdp.s(1,1),(12 * (j-1) +1): (12*j))
            text(text_coords(j,1),text_coords(j,2),choice{j},'FontSize',14,'FontWeight','Bold','color','blue')
        else
            text(text_coords(j,1),text_coords(j,2),choice{j},'FontSize',14)
        end
    end
    
    hold on;
    
    if t_high_i < T_Higher
        
        for t_sacc = quadrant_arrival_times(t_high_i):(quadrant_arrival_times(t_high_i+1))
            
            if t_sacc+1 > size(T,1)
                break
            else
                plot(T(t_sacc:t_sacc+1,1),T(t_sacc:t_sacc+1,2),'b--','LineWidth',3); hold on;
                pause(0.05);
            end
            
            M(frame_num) = getframe(gcf);
            frame_num = frame_num + 1;
            pause(0.05);            
            
        end
        
        if mdp.s(2,t_high_i+1) > 5
            
            subplot(121)
            
            switch mdp.o(2,t_high_i+1)
                case 1
                    text(-2,-11,'UNDECIDED','FontSize',20,'FontWeight','Bold','color',[0.5 0.5 0.5])
                case 2
                    text(-2,-11,'CORRECT!','FontSize',20,'FontWeight','Bold','color',[0.1 0.8 0.2])
                case 3
                    text(-2,-11,'INCORRECT!','FontSize',20,'FontWeight','Bold','color',[0.8 0.1 0.4])
            end
            
            subplot(122)
            
            switch mdp.o(2,t_high_i+1)
                case 1
                    text(-2,-11,'UNDECIDED','FontSize',20,'FontWeight','Bold','color',[0.5 0.5 0.5])
                case 2
                    text(-2,-11,'CORRECT!','FontSize',20,'FontWeight','Bold','color',[0.1 0.8 0.2])
                case 3
                    text(-2,-11,'INCORRECT!','FontSize',20,'FontWeight','Bold','color',[0.8 0.1 0.4])
            end
            
            M(frame_num) = getframe(gcf);
            frame_num = frame_num + 1;
            pause(0.05); 
            
        end

  
    %% parameterize your dots 
    
        if mdp.mdp(t_high_i+1).s(1,1) > 1 % only spend resources drawing dots this if there's actually a dot pattern at the next quadrant

            % center around which to seed random dots

            quad_idx = mdp.o(3,t_high_i+1)-1;

            center = quadrant_centers(quad_idx,:);

            % seed dots around center and with a spread proportional to apertureSize
            x = (rand(1,nDots)-.5) * apertureSize(1) + center(1);
            y = (rand(1,nDots)-.5) * apertureSize(2) + center(2);

            switch mdp.mdp(t_high_i+1).s(1,1)
                case 2
                    direction = 180;
                case 3
                    direction = 90;
                case 4
                    direction = 0;
                case 5
                    direction = 270;
            end

            coherence = mdp.mdp(t_high_i+1).A{1}(2,2,1);

            % get indices of dots that will be moving coherently
            coh_idx = rand(1,nDots) < coherence;

            % initialize motion vectors of all dots to random direction
            rand_directions = 90.*randi(4,1,nDots) .* pi/180;
%             rand_directions = 360 .* rand(1,nDots) .* pi/180;

            % change random directions into actual increments in x/y position
            dx = 0.2*sin(rand_directions);
            dy = -0.2*cos(rand_directions);

            % change coherent directions of coh_idx dots into actual increments in
            % x/y position
            dx(coh_idx) = 0.2*sin(direction * pi/180);
            dy(coh_idx) = -0.2*cos(direction * pi/180);

            l = center(1)-apertureSize(1)/2;
            r = center(1)+apertureSize(1)/2;
            b = center(2)-apertureSize(2)/2;
            t = center(2)+apertureSize(2)/2;

            lives = ceil(rand(1,nDots).*lifetime);

            sacc_duration = mdp.mdp(t_high_i+1).T;

            numSamples = sacc_duration * 5;

            dotHandles = [];

            for samp = 1:numSamples

                subplot(122)

                if samp > 1
                    delete(dotHandles(samp-1));
                end

                x = x + dx;
                y = y + dy;

                x(x<l) = x(x<l) + apertureSize(1);
                x(x>r) = x(x>r) - apertureSize(1);
                y(y<b) = y(y<b) + apertureSize(2);
                y(y>t) = y(y>t) - apertureSize(2);

                % increment the lives of the dots
                lives = lives + 1;

                %find the 'dead' dots
                deadDots = mod(lives,lifetime)==0;

                %replace the positions of the dead dots to a random location
                x(deadDots) = (rand(1,sum(deadDots))-.5).*apertureSize(1) + center(1);
                y(deadDots) = (rand(1,sum(deadDots))-.5).*apertureSize(2) + center(2);

                hold on;
                dotHandles(samp) = scatter(x,y,'ko','MarkerFaceColor','w');

                if mod(samp,5) == 0

                    LL_s = zeros(size(UP)); % posterior estimate of the lower-level hidden state

                    for o_i = 1:5
                        p = mdp.mdp(t_high_i+1).xn{1}(end,o_i,samp/5,samp/5);
                        LL_s = LL_s + p.*eval(outcomes{o_i});
                    end

                    LL_s = LL_s./max(LL_s(:));

                    subplot(121)

                    quad_idx = mdp.o(3,t_high_i+1) - 1;

                    imagesc([center(1) - 2.5, center(1) + 2.5], ...
                        [center(2) - 2.5, center(2) + 2.5],flipud(LL_s)); % have to flip because of call 'axis image xy' (line 137 above)


                end

                pause(0.01);

                M(frame_num) = getframe(gcf);
                dot_idx = [dot_idx,frame_num];
                frame_num = frame_num + 1;
                

            end

            % clear the dots at the end

            delete(dotHandles(end));
            
        end
    
    end
    
   
    M(frame_num) = getframe(gcf);
    frame_num = frame_num + 1;
    
    pause(0.1);
    
    
end
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
end












