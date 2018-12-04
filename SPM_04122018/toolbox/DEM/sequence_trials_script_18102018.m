clear all; clc;

[ all_configs,scene_idx ] = generate_scenes();

% uncertainty = 20; % just blow it up so it's pretty uncertain what the dot pattern is
noise_vals = linspace(0,2,16); 

prior_scene_probs = linspace(0.25,0.8,16);

N = 10;

% all_results = cell(1,length(prior_scene_probs));
all_results = cell(length(noise_vals),length(prior_scene_probs));

save('MDPresults_EV_24102018.mat','all_results','noise_vals','prior_scene_probs')

%% run it

for i = 1:length(noise_vals)
    
    for j = 1:length(prior_scene_probs)
    
    
        % randomly initialize which scene the agent believes in
        prior_scene_beliefs = cell(N,1);
        idx_tmp = randi([1 4],N,1);
        prior_scene_beliefs(idx_tmp == 1) = {'UP_RIGHT'};
        prior_scene_beliefs(idx_tmp == 2) = {'RIGHT_DOWN'};
        prior_scene_beliefs(idx_tmp == 3) = {'DOWN_LEFT'};
        prior_scene_beliefs(idx_tmp == 4) = {'LEFT_UP'};
        clear idx_tmp;
        
        % randomly initialize true hidden states
        true_scenes = randi([1 48],N,1);
        
        for trial = 1:N
            
            MDPdeep = initialize_MDPdeep(prior_scene_beliefs{trial},prior_scene_probs(j),all_configs,scene_idx,true_scenes(trial));
            MDPshallow = initialize_MDPshallow_witha(noise_vals(i));
            % nest the shallow MDP within the deeper one
            MDPdeep.MDP = MDPshallow;
            MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));
            
            MDPfull(1,trial) = spm_MDP_check(MDPdeep);
            clear MDPdeep MDPshallow;
            
        end
        
        %% Solving
        
        tic
        % solve a sequence of trials
        %==========================================================================
        MDPresult  = spm_MDP_VB_X(MDPfull); clear MDPfull;
        all_results{i,j} = MDPresult;
        
        clear MDPresult;
        
        fprintf('Time taken to complete %d trials; %.2f minutes\n',N,toc/60)
        
    end
end


%% get the epistemic value over time from the quadrants

load('/Users/conorheins/Documents/MATLAB/spm12_r6906_mexmaci64/toolbox/DEM/MDPresults_EV_24102018.mat')

%%
[all_configs, scene_idx] = generate_scenes();
var_names = {'EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','Time to Choose','Uncertainty','Prior precision','Accuracy','CorrectPrior'};

stats_matrix = zeros(16*16*10,length(var_names));

row_idx = 1;
for i = 1:length(noise_vals)
    
    for j = 1:length(prior_scene_probs)
        
        for n_i = 1:10
            
            current_trial = all_results{i,j}(n_i);
            stats_matrix(row_idx,1:8) = mean(current_trial.EV(2:5,:),1);
            time2choose = find(current_trial.o(3,:) > 5,1)-1;
            if isempty(time2choose)
                stats_matrix(row_idx,9) = size(current_trial.o,2);
            else
                stats_matrix(row_idx,9) = time2choose;
            end
            stats_matrix(row_idx,10) = noise_vals(i);
            stats_matrix(row_idx,11) = prior_scene_probs(j);
            
            % second row of o has the accuracy read-out: 1 if correct, 2 if
            % incorrect, 3 if undecided
            o2 = current_trial.o(2,:);
            if any(o2 == 2) && ~any(o2 == 3)
                stats_matrix(row_idx,12) = 1;
            elseif ~any(o2 == 2 | o2 == 3)
                stats_matrix(row_idx,12) = 3;
            else
                stats_matrix(row_idx,12) = 2;
            end
            
            true_state = current_trial.s(1,1);
            believed_state = find(current_trial.D{1},1);
            
            for scene_i = 1:4
                if and(ismember(true_state,scene_idx{scene_i}),ismember(believed_state,scene_idx{scene_i}))
                    stats_matrix(row_idx,13) = 1;
                    break
                else
                    stats_matrix(row_idx,13) = 0;
                end
            end
            
            row_idx = row_idx + 1;
        end
    end
end

% do some plots bro

% plots of average epistemic value over time, for trials last 5,6,7,8
% saccades long, separated by noise values

i = 1;
for trial_length = 5:8

    figure(i);
    idx = ismember(stats_matrix(:,10),noise_vals(1:4));
    plot(mean(stats_matrix(idx&stats_matrix(:,9) == trial_length,1:trial_length),1),'DisplayName','Low Noise')
    
    idx = ismember(stats_matrix(:,10),noise_vals(5:8));
    hold on; plot(mean(stats_matrix(idx&stats_matrix(:,9) == trial_length,1:trial_length),1),'DisplayName','MedLow Noise')
    
    idx = ismember(stats_matrix(:,10),noise_vals(9:12));
    hold on; plot(mean(stats_matrix(idx&stats_matrix(:,9) == trial_length,1:trial_length),1),'DisplayName','MedHigh Noise')
    
    idx = ismember(stats_matrix(:,10),noise_vals(13:16));
    hold on; plot(mean(stats_matrix(idx&stats_matrix(:,9) == trial_length,1:trial_length),1),'DisplayName','High Noise')
    title(sprintf('Epistemic value and noise for %d-timestep trials',trial_length))
    legend('show')
    
    i = i + 1;

end

% plots of average epistemic value over time, for trials last 5,6,7,8
% saccades long, separated by prior precisions

i = 1;
for trial_length = 5:8

    figure(i);
    idx = and(ismember(stats_matrix(:,11),prior_scene_probs(1:4)),stats_matrix(:,13) == 1);
    plot(mean(stats_matrix(idx&stats_matrix(:,9) == trial_length,1:trial_length),1),'DisplayName','Low Precision')
    
    idx = and(ismember(stats_matrix(:,11),prior_scene_probs(5:8)),stats_matrix(:,13) == 1);
    hold on; plot(mean(stats_matrix(idx&stats_matrix(:,9) == trial_length,1:trial_length),1),'DisplayName','MedLow Precision')
    
    idx = and(ismember(stats_matrix(:,11),prior_scene_probs(9:12)),stats_matrix(:,13) == 1);
    hold on; plot(mean(stats_matrix(idx&stats_matrix(:,9) == trial_length,1:trial_length),1),'DisplayName','MedHigh Precision')
    
    idx = and(ismember(stats_matrix(:,11),prior_scene_probs(13:16)),stats_matrix(:,13) == 1);
    hold on; plot(mean(stats_matrix(idx&stats_matrix(:,9) == trial_length,1:trial_length),1),'DisplayName','High Precision')
    title(sprintf('Epistemic value and prior precision for %d-timestep trials',trial_length))
    legend('show')
    
    i = i + 1;

end

%%

% find all trials with same hidden state configurations but different prior beliefs 
% and compare evolution of beliefs/epistemic value for them


        
        
        
        
        