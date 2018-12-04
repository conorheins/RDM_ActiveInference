function [ M ] = MDP_beliefs_video_ENIPC(MDP,trial_idx,all_configs)
% MDP_BELIEFS_VIDEO Creates video of the trajectory of an MDP active inference agent's
% visual search behavior during a single trial, as well as the
% corresponding evolution of posterior beliefs (at multiple hierarchical
% levels)

mdp = MDP(1,trial_idx);

outcomes = {'null','UP','RIGHT','DOWN','LEFT'};

scene_config    = all_configs(mdp.s(1,1),:);
scene_config(strcmp(scene_config,'n')) = outcomes(1);
scene_config(strcmp(scene_config,'U')) = outcomes(2);
scene_config(strcmp(scene_config,'R')) = outcomes(3);
scene_config(strcmp(scene_config,'D')) = outcomes(4);
scene_config(strcmp(scene_config,'L')) = outcomes(5);

figure('Position',[300 300 800 600]);

% load images
%----------------------------------------------------------------------
load MDP_search_graphics_DOT
null = zeros(size(UP)) + 1;

%--------------------------------------------------------------------------
Ni    = 1:size(mdp.xn{1},1);            % index vector for variational updates
Nx    = length(Ni);                     % total number of variational updates
Ne = find([mdp.o(3,:) > 5,1],1);        % number of eye-movements until categorization/scene choice (including scene choice itself)


% first subplot shows true scene & entailed outcomes
subplot(1,2,1)

% locations for first subplot (includes the fixation center, the four quadrants
% and the four category choices
%--------------------------------------------------------------------------
choice_spacing = linspace(-1.25,1.1,4);
x1 = [0,0;-1 -1; -1 1; 1 -1;1 1;-1.5 2.5;-0.75 2.5;0.2,2.5;1.2 2.5];
x1(6:end,1) = choice_spacing;
y  = x1 + 1/6;
r1 = [-1,1]/1.5;

displace_idx = find(~strcmp(scene_config,'null')) + 1;
y(displace_idx,:) = y(displace_idx,:) + 1/5;

for i = 1:numel(scene_config)
    image(r1 + x1(i + 1,1),r1 + x1(i + 1,2),eval(scene_config{i})), hold on
end

% choices
%----------------------------------------------------------------------
choice = {'UR','RD','DL','LU'};
for i = 1:4
    if ismember(mdp.s(1,1),(12 * (i-1) +1): (12*i))
        text(x1(i+5,1),x1(i+5,2),choice{i},'FontSize',14,'FontWeight','Bold','color','red')
    else
        text(x1(i+5,1),x1(i+5,2),choice{i},'FontSize',14)
    end
end

% labels
%----------------------------------------------------------------------
for i = 1:size(x1,1)
    if ismember(i,[6,7,8,9])
        text(y(i,1)-0.1,y(i,2),num2str(i),'FontSize',16,'FontWeight','Bold','color',[0.2 0.5 0.8])
    else
        text(y(i,1)-0.25,y(i,2),num2str(i),'FontSize',20,'FontWeight','Bold','color',[0.2 0.8 0.8])
    end
end

axis([-2,2,-2,3])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

for i = 1:numel(mdp.o(3,:))
    X(i,:) = x1(mdp.o(3,i),:);
end
for j = 1:2
    T(:,j) = interp(X(:,j),8,2);
    T(:,j) = T(:,j) + spm_conv(0.8*randn(size(T(:,j))),2)/16;
end

quadrant_arrival_times = [1,8:8:length(T)];


% second subplot shows fixation center & beliefs about scene (projected
% onto their expected on the four quadrants )

% locations for the second subplot (includes a central location indicating
% fixational outcome, and four quadrants indicating beliefs about outcomes
% of states at each of those four quadrants)
%--------------------------------------------------------------------------
x2 = x1(1:5,:);
r2 = [-1 1]/2;
UP = double(imresize(UP,0.8));
RIGHT = double(imresize(RIGHT,0.8));
DOWN = double(imresize(DOWN,0.8));
LEFT = double(imresize(LEFT,0.8));

% make a mask that represents currently-foveated outcome
null = zeros(size(UP)) + 1;
mask  = hamming(size(UP,1));
mask  = mask*mask';
for i = 1:3
    mask(:,:,i) = mask(:,:,1);
end

frame_num = 1;
for sacc_i = 1:Ne
    
    subplot(1,2,1)
    title(sprintf('Timestep # %d',sacc_i))
    
    subplot(1,2,2)
    title(sprintf('Timestep # %d',sacc_i))
    axis([-2,2,-2,2])
    axis image ij
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    if sacc_i > size(X,1)
        subplot(122)
        text(-.5,0,'UNDECIDED','FontSize',20,'FontWeight','Bold','color',[0.5 0.5 0.5])
        M(frame_num) = getframe(gcf);
        frame_num = frame_num + 1;         
        break
    end
    
    subplot(1,2,1)
    plot(X(sacc_i,1),X(sacc_i,2),'bo','MarkerSize',20,'LineWidth',2);
    
    if sacc_i > Ne
        plot(T(1:end,1),T(1:end,2),'r:','LineWidth',3);
        M(frame_num) = getframe(gcf);
        frame_num = frame_num + 1;
    elseif sacc_i == 1
        plot(T(1,1),T(1,2),'r:','LineWidth',3);
        M(frame_num) = getframe(gcf);
        frame_num = frame_num + 1;
    else
        for t = quadrant_arrival_times(sacc_i-1):quadrant_arrival_times(sacc_i)
            plot(T(t:t+1,1),T(t:t+1,2),'r:','LineWidth',3); hold on;
            M(frame_num) = getframe(gcf);
            frame_num = frame_num + 1;
            pause(0.05);
        end
    end
    
    
    subplot(1,2,2)
    

    hold on;
    
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
        p     = mdp.xn{1}(end,config_i,sacc_i,sacc_i);
        
        for j = 1:4
            S{j} = S{j} + p.*eval(a{j});
        end
        
    end
    
    if exist('h','var')
        delete(h);
    end
    for j = 1:numel(S)
        h(j) = imagesc(r2 + x2(j + 1,1),r2 + x2(j + 1,2),S{j}/max(S{j}(:)));
    end
    
    d   = (1 - exp(1 - Nx));
    
    LO = zeros(size(UP));
 
    for o_i = 1:5
        LO = LO + mdp.O{1,sacc_i}(o_i)*eval(outcomes{o_i});
    end
    
    LO = LO./max(LO(:));
    
    imagesc(r2,r2,LO.*mask*d);
    
    if sacc_i == Ne
        switch mdp.o(2,sacc_i)
            case 1
                text(-.5,0,'UNDECIDED','FontSize',20,'FontWeight','Bold','color',[0.5 0.5 0.5])
            case 2
                text(-.5,0,'CORRECT!','FontSize',20,'FontWeight','Bold','color',[0.1 0.8 0.2])
            case 3
                text(-.5,0,'INCORRECT!','FontSize',20,'FontWeight','Bold','color',[0.8 0.1 0.4])
        end
    end
    
    title(sprintf('Timestep # %d',sacc_i))
    M(frame_num) = getframe(gcf);
    frame_num = frame_num + 1;
    
    pause(1.5);
    
end
    


end

