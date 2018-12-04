function view_trial_trajectory(MDP_structure,trial_idx,Scenes)
%VIEW_TRIAL_TRAJECTORY Given an MDP structure, chooses a trial (given by
% trial_idx) and displays the saccade trajectory superimposed on the visual
% display 

% in the absence of 'trial_idx', random trial is selected
if ~exist('trial_idx','var') || isempty(trial_idx)
    trial_to_plot  = MDP_structure(randi(length(MDP_structure),1));
else
    trial_to_plot = MDP_structure(trial_idx);
end

% display status of the visual display, given hidden state configuration
Scene_config = Scenes(trial_to_plot.s(1,1),:);
displaySceneConfig(Scene_config,[],0.04);
    
%draw saccades onto existing display

hold on;
xlimz = xlim;
ylimz = ylim;
for sacc_i = 1:size(trial_to_plot.o,2)
    
    switch trial_to_plot.o(2,sacc_i)
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
        case 7 % choose scene B (scene B 'button')
            center_coords = [7*(xlimz(2)/16), ylimz(2)/8];
        case 8 % choose scene C (scene C 'button')
            center_coords = [9*(xlimz(2)/16), ylimz(2)/8];
        case 9 % choose scene D (scene D 'button')
            center_coords = [11*(xlimz(2)/16), ylimz(2)/8];
    end
    
    if exist('last_coords','var')
        plot([last_coords(1) center_coords(1)],[last_coords(2) center_coords(2)],'r--','LineWidth',5)
    end
    plot(center_coords(1),center_coords(2),'ro','LineWidth',10,'MarkerSize',30);
    text(center_coords(1)+10,center_coords(2)+10,num2str(sacc_i),'FontSize',26,'FontWeight','bold')
    
    last_coords = center_coords;
    
    if trial_to_plot.o(2,sacc_i) >= 6
        if trial_to_plot.o(1,sacc_i) == 6
            text(center_coords(1),center_coords(2) - 20, 'CORRECT!','FontSize',30,'Color',[0 0.5 0.5])
        elseif trial_to_plot.o(1,sacc_i) == 7
            text(center_coords(1),center_coords(2) - 20, 'INCORRECT!','FontSize',30,'Color',[0.5 0 0.1])
        end
        break
    end
end



end