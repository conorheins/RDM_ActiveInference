%%
clear all; clc;

[ all_configs,scene_idx ] = generate_scenes();

noise_vals = exp(linspace(log(1),log(5),20));

prior_scene_probs = linspace(0.25,0.75,16);

N = 250; % number of trials for each {noise_val(i), prior_scene_probs(j)} combination
T = 32; % temporal depth of lower level

all_results = cell(length(noise_vals),length(prior_scene_probs));

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
            MDPshallow = initialize_MDPshallow_witha_newp(noise_vals(i),T);
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

%% analyze

stats_matrix= zeros([size(all_results),2]);

for i = 1:length(noise_vals)
    
    for j = 1:length(prior_scene_probs)
        
        
        trial_stats = zeros(1,N);
        
        for trial = 1:N
            
            current_data = all_results{i,j}(trial);
            
            saccades2inspect = find(current_data.o(1,:) > 1);
            
            for s_i = 1:length(saccades2inspect)
                trial_stats(trial) = trial_stats(trial)  + size(current_data.mdp(saccades2inspect(s_i)).X{1},2)/length(saccades2inspect);
            end
            
        end
       
        stats_matrix(i,j,1) = mean(trial_stats);
        stats_matrix(i,j,2) = std(trial_stats);
       
        
    end
end



