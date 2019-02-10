function [ MDP ] = random_searchMDP(prior_belief,prior_precis,true_scene,sens_precis)
%RANDOM_SEARCHMDP Runs a single instance of hierarchical scene construction MDPs
%and spits them out into MDP result

T = 20; % number of timesteps at the lower level
policy_depth = 1; % depth of policies at the lower level
chi_val = 1/10e30; % bound on posterior negentropy of beliefs at the lower level 

[ all_configs,scene_idx ] = generate_scenes();

MDPdeep = initialize_MDPdeep(prior_belief,prior_precis,all_configs,scene_idx,true_scene);
MDPshallow = initialize_MDPshallow_v7(sens_precis,T,policy_depth);
MDPshallow.chi = chi_val;

% nest the shallow MDP within the deeper one
MDPdeep.MDP = MDPshallow;
MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));

MDPfull = spm_MDP_check(MDPdeep);
clear MDPdeep MDPshallow;

%% Solving

% solve a sequence of trials
%==========================================================================
MDP  = spm_MDP_VB_X_RCH(MDPfull);

end

