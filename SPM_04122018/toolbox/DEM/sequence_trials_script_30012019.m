%%

T = 20; % number of timesteps at the lower level
policy_depth = 1; % depth of policies at the lower level
chi_val = 1/10e30; % bound on posterior negentropy of beliefs at the lower level 

results_array = []; % initialize empty array to contain results

iter = 1;
while iter <= 80
% while iter <= 2
    %% solve and evaluate performance for sequence of trials
    
    tic

    [ all_configs,scene_idx ] = generate_scenes();
    
    N = 300; % number of trials
%     N = 3; % debugging

    % randomly initialize agents' beliefs
    prior_scene_beliefs = cell(N,1);
    idx_tmp = randi([1 4],N,1);
    prior_scene_beliefs(idx_tmp == 1) = {'UP_RIGHT'};
    prior_scene_beliefs(idx_tmp == 2) = {'RIGHT_DOWN'};
    prior_scene_beliefs(idx_tmp == 3) = {'DOWN_LEFT'};
    prior_scene_beliefs(idx_tmp == 4) = {'LEFT_UP'};
    
    % randomly initialize agents' precision in their beliefs
    prior_scene_probs = 0.25 + (0.8 - 0.25).*rand(N,1);
    
    % randomly initialize true hidden states
    true_scenes = randi([1 48],N,1);
    
    % randomly initialize sensory uncertainty
    noise_vals = 5.5.*rand(N,1); 
    
    for trial = 1:N
        
        MDPdeep = initialize_MDPdeep(prior_scene_beliefs{trial},prior_scene_probs(trial),all_configs,scene_idx,true_scenes(trial));
        MDPshallow = initialize_MDPshallow_v7(noise_vals(trial),T,policy_depth);
        MDPshallow.chi = chi_val;
        
        % nest the shallow MDP within the deeper one
        MDPdeep.MDP = MDPshallow;
        MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));
        
        MDPfull(1,trial) = spm_MDP_check(MDPdeep);
        clear MDPdeep MDPshallow;
        
    end
    
    %% Solving
    
    % solve a sequence of trials
    %==========================================================================
    MDPresult  = spm_MDP_VB_X_RCH(MDPfull); clear MDPfull;
    
    % gather trial-wise info (time-course of H_int as well
                        
    for trial = 1:N
        
        curr_trial = MDPresult(trial);
        
        quadrant_saccades = find(and(curr_trial.o(1,:) > 1,curr_trial.o(1,:) < 6));
        
        trial_array = [];
        
        for ii = 1:length(quadrant_saccades)
            
            sacc_ii = curr_trial.mdp(quadrant_saccades(ii));
            
            sub_array = zeros(1,2); % one extra column (on the left) to store saccade number
            
            sub_array(1,1) = quadrant_saccades(ii);
            
            break_time = size(sacc_ii.o,2); %the first time they made a decision at the lower level ('perceptual policy')
            sub_array(1,2) = break_time;

            trial_array = [trial_array; [trial, sub_array]];
            
        end
        
        results_array = [results_array;[noise_vals(trial)*ones(size(trial_array,1),1),prior_scene_probs(trial)*ones(size(trial_array,1),1),...
            true_scenes(trial)*ones(size(trial_array,1),1), idx_tmp(trial)*ones(size(trial_array,1),1), trial_array]];
        
    end
    
    save(sprintf('W:/#Common/4Conor/RDPsearch_Hierarch15_results/RDPsearch_Hierarch15_%d_%s.mat',iter,datestr(now,'ddmmyy')),'chi_val','policy_depth','T',...
       'MDPresult','all_configs','scene_idx','prior_scene_beliefs','prior_scene_probs','true_scenes','noise_vals','-v7.3')

    clear MDPresult;
    
    fprintf('Time taken to complete iteration %d (%d trials); %.2f minutes\n',iter,N,toc/60)
    
    iter = iter + 1;

    
end

save(sprintf('W:/#Common/4Conor/RDPsearch_Hierarch15_results/RDPsearch_Hierarch15_%s_results.mat',datestr(now,'ddmmyy')),'chi_val','policy_depth','T',...
       'results_array','all_configs','scene_idx','-v7.3');

clear all; 


