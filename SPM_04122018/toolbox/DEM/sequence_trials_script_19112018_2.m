%%
iter = 1;
while iter <= 80
% while iter <= 1
    %% solve and evaluate performance for sequence of trials
    
    tic

    [ all_configs,scene_idx ] = generate_scenes();
    
    N = 300; % number of trials
%     N = 10; % debugging
    
    % randomly initialize agents' beliefs
    prior_scene_beliefs = cell(N,1);
    idx_tmp = randi([1 4],N,1);
    prior_scene_beliefs(idx_tmp == 1) = {'UP_RIGHT'};
    prior_scene_beliefs(idx_tmp == 2) = {'RIGHT_DOWN'};
    prior_scene_beliefs(idx_tmp == 3) = {'DOWN_LEFT'};
    prior_scene_beliefs(idx_tmp == 4) = {'LEFT_UP'};
    clear idx_temp;
    
    % randomly initialize agents' precision in their beliefs
    prior_scene_probs = 0.25 + (0.8 - 0.25).*rand(N,1);
    
    % randomly initialize true hidden states
    true_scenes = randi([1 48],N,1);
    
    % randomly initialize sensory uncertainty
    noise_vals = 5.5.*rand(N,1); 
    
    for trial = 1:N
        
        MDPdeep = initialize_MDPdeep(prior_scene_beliefs{trial},prior_scene_probs(trial),all_configs,scene_idx,true_scenes(trial));
        MDPshallow = initialize_MDPshallow_witha_newp(noise_vals(trial),32);
        % nest the shallow MDP within the deeper one
        MDPdeep.MDP = MDPshallow;
        MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));
        
        MDPfull(1,trial) = spm_MDP_check(MDPdeep);
        clear MDPdeep MDPshallow;
        
    end
    
    %% Solving
    
    % solve a sequence of trials
    %==========================================================================
    MDPresult  = spm_MDP_VB_X(MDPfull); clear MDPfull;
    
   save(sprintf('W:/#Common/4Conor/RDPsearch_Hierarch8_results/RDPsearch_Hierarch8_%d_%s.mat',iter,datestr(now,'ddmmyy')),'MDPresult','all_configs',...
        'scene_idx','prior_scene_beliefs','prior_scene_probs','true_scenes','noise_vals','-v7.3')

   
    clear MDPresult;
    
    fprintf('Time taken to complete iteration %d (%d trials); %.2f minutes\n',iter,N,toc/60)
    
    iter = iter + 1;

    
end

clear all; 


