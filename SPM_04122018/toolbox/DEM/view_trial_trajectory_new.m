function view_trial_trajectory_new(MDP,trial_idx,Scenes,video_mode)

if ~exist('video_mode','var') || isempty(video_mode)
    video_mode = false;
end

mdp = MDP(1,trial_idx);

% locations
%--------------------------------------------------------------------------
x = [0,0;-1 -1; -1 1; 1 -1;1 1;-1.25,2.5;-0.4167 ,2.5;0.4167,2.5;1.25,2.5];
y = x + 1/6;
r = [-1,1]/1.25;

% plot cues
%--------------------------------------------------------------------------

figure('Position',[400 400 550 700]);
% load images
%----------------------------------------------------------------------
load MDP_search_graphics_DOT
null = zeros(size(UP)) + 1;
scene_config    = Scenes(mdp.s(1,1),:);
scene_config(strcmp(scene_config,'U')) = {'UP'};
scene_config(strcmp(scene_config,'R')) = {'RIGHT'};
scene_config(strcmp(scene_config,'D')) = {'DOWN'};
scene_config(strcmp(scene_config,'L')) = {'LEFT'};
scene_config(strcmp(scene_config,'n')) = {'null'};

displace_idx = find(~strcmp(scene_config,'null')) + 1;
y(displace_idx,:) = y(displace_idx,:) + 1/5;


for i = 1:numel(scene_config)
    image(r + x(i + 1,1),r + x(i + 1,2),eval(scene_config{i})), hold on
end

% choices
%----------------------------------------------------------------------
choice = {'UP RIGHT','RIGHT DOWN','DOWN LEFT','LEFT UP'};
for i = 1:4
    if ismember(mdp.s(1,1),(12 * (i-1) +1): (12*i))
        text(y(i+5,1)-0.4,2.5,choice{i},'FontSize',12,'FontWeight','Bold','color','red')
    else
        text(y(i+5,1)-0.4,2.5,choice{i},'FontSize',12)
    end
end
axis image, axis([-2,2,-2,3])

% labels
%----------------------------------------------------------------------
for i = 1:size(x,1)
    if ismember(i,[6,7,8,9])
        text(y(i,1)-0.15,y(i,2),num2str(i),'FontSize',16,'FontWeight','Bold','color',[0.2 0.5 0.8])
    elseif i == 1
        text(y(i,1)-0.25,y(i,2),num2str(i),'FontSize',20,'FontWeight','Bold','color',[0.2 0.8 0.8])
    else
        text(y(i,1),y(i,2),num2str(i),'FontSize',20,'FontWeight','Bold','color',[0.2 0.8 0.8])
    end
end


% Extract and plot eye movements
%--------------------------------------------------------------------------
for i = 1:numel(mdp.o(3,:))
    X(i,:) = x(mdp.o(3,i),:);
end
for j = 1:2
    T(:,j) = interp(X(:,j),8,2);
    T(:,j) = T(:,j) + spm_conv(0.8*randn(size(T(:,j))),2)/16;
end
if video_mode
    quadrant_arrival_times = [1,8:8:length(T)];
    saccade_num = 1;
    for t = 1:length(T)
        if ismember(t,quadrant_arrival_times)
            if saccade_num <= size(X,1)
                plot(X(saccade_num,1),X(saccade_num,2),'bo','MarkerSize',20,'LineWidth',2);
                saccade_num = saccade_num + 1;
            end
        end
        plot(T(1:t,1),T(1:t,2),'r:','LineWidth',3);
        pause(0.1);
    end
else
    plot(T(:,1),T(:,2),'r:','LineWidth',3)
    plot(X(:,1),X(:,2),'bo','MarkerSize',20,'LineWidth',2)
end
