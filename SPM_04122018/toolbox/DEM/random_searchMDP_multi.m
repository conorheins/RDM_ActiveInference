function [ MDP ] = random_searchMDP_multi(N,max_prior_b,max_precis)
%RANDOM_SEARCHMDP Just runs nTrials of hierarchical scene construction MDPs
%and spits them out into MDP result
% INPUTS: nTrials -- how many trials in the sequence

if nargin < 3
    max_precis = 1;
end

if nargin < 2
    max_prior_b = 0.7;
end

if nargin < 1
    N = 1;
end

T = 20; % number of timesteps at the lower level
policy_depth = 1; % depth of policies at the lower level
chi_val = 1/10e30; % bound on posterior negentropy of beliefs at the lower level 

[ all_configs,scene_idx ] = generate_scenes();

% randomly initialize agents' beliefs
prior_scene_beliefs = cell(N,1);
idx_tmp = randi([1 4],N,1);
prior_scene_beliefs(idx_tmp == 1) = {'UP_RIGHT'};
prior_scene_beliefs(idx_tmp == 2) = {'RIGHT_DOWN'};
prior_scene_beliefs(idx_tmp == 3) = {'DOWN_LEFT'};
prior_scene_beliefs(idx_tmp == 4) = {'LEFT_UP'};

% randomly initialize agents' precision in their beliefs
prior_scene_probs = 0.25 + (max_prior_b - 0.25).*rand(N,1);

% randomly initialize true hidden states
true_scenes = randi([1 48],N,1);

% randomly initialize sensory uncertainty
noise_vals = max_precis.*rand(N,1);

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
MDP  = spm_MDP_VB_X_RCH(MDPfull);

end

