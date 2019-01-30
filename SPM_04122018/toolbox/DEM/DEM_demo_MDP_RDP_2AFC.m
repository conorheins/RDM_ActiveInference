%% DEM_demo_MDP_RDPsearch
%  R. Conor Heins
%
%   MDP formulation of a two alternative forced choice task  with
%   2 different random dot patterns (hidden states), moving in orthogonal 
%   directions (e.g. LEFT and RIGHT).
%
%   A hierarchical MDP formulation is used in order to enable 
%   introduction of sensory uncertainty at the lower hierarchical level.

% 1. The first (shallowest) hierarchical level describes the mapping from
%    a 'true' individual RDP-direction (e.g. a LEFT RDP) to the observed/sampled direction. 
%    In the trivial case of no sensory uncertainty, this mapping is a bijective 
%    identity mapping, so that an RDP-direction at the higher level is directly 
%    observed as such at the lower level. The likelihood mapping at this first level,
%    from true-direction to observed direction, can however be manipulated as a model
%    of sensory uncertainty e.g. incoherence of dot-motion in RDPs.
% 2. The second (deepest) hierarchical level describes the mapping from the latent scene 
%    (e.g. 'LEFT') to the dot pattern at the middle hierarchical level (the
%    hidden scenes for the first level.

% $Id: DEM_demo_MDP_RDP_2AFC.m  2019-01-20 Conor Heins $

clear all; close all;
%% 1. initialize the first level: 'shallowest' MDP -- mapping RDP directions to (optionally noisy versions of) themselves in a given quadrant

% noise parameter, manipulates the coherence of dot motion via the inverse temperature parameter of a softmax transformation of
% a categorical distribution over dot directions, with a 1 at the true dot motion direction
p = 0.25; 

% temporal depth of updates at the lower level
T = 32; 

% depth of 'policies' at lower level
policy_depth = 1;

MDPshallow = initialize_MDPshallow_2AFC(p,T,policy_depth);

%% 2. initialize the second level: the'deepest' MDP 
%     maps hidden state (LEFT/RIGHT) to RDP manifestation (also simply LEFT / RIGHT) in the single FOV

% choose true scene identity 
true_RDP = 'LEFT';

% initialize highest level of the MDP
MDPdeep = initialize_MDPdeep_2AFC(true_RDP);

%% 3. Link the two levels into a single hierarchical model

% nest the shallow MDP within the deeper one
MDPdeep.MDP = MDPshallow; 
MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));

MDPfull = spm_MDP_check(MDPdeep);

%% 4. Solving

% illustrate a single trial
%==========================================================================
MDPresult = spm_MDP_VB_X(MDPfull);

%% 5. Run a bunch of trials with different sensory precisions

% temporal depth of updates at the lower level
T = 32; 

% depth of 'policies' at lower level
policy_depth = 1;

% choose true scene identity 
true_RDP = 'LEFT';

P_vec = [0.25 0.5 0.8 1 1.2 2];

N = 50;

xn_array = [];
H_int_array = [];
rt_array = [];

for p_i = 1:length(P_vec)
    
    
    MDPshallow = initialize_MDPshallow_2AFC(P_vec(p_i),T,policy_depth);
    MDPdeep = initialize_MDPdeep_2AFC(true_RDP);
    MDPdeep.MDP = MDPshallow;
    MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));
    
    for trial = 1:N        
        MDPfull(1,trial) = spm_MDP_check(MDPdeep);
    end
    
    clear MDPdeep MDPshallow;
    
    MDPresult = spm_MDP_VB_X(MDPfull);
    
    for trial = 1:N
        
        xn = squeeze(MDPresult(trial).mdp(2).xn{1}(end,:,:,:));
        
        both_beliefs = zeros(1,2,32);
        
        both_beliefs(1,1,1:size(xn,3))  = diag(squeeze(xn(2,:,:)))';
        both_beliefs(1,2,1:size(xn,3))  = diag(squeeze(xn(3,:,:)))';
        
        xn_array = [xn_array ; both_beliefs];
        
        temp_Hint_vec = zeros(1,32);
        temp_Hint_vec(1:length(MDPresult(trial).mdp(2).H_int)) = MDPresult(trial).mdp(2).H_int;
        H_int_array = [H_int_array;[P_vec(p_i),temp_Hint_vec]];
        
        temp_rt_vec = zeros(1,32);
        temp_rt_vec(1:length(MDPresult(trial).mdp(2).rt)) = MDPresult(trial).mdp(2).rt;
        rt_array = [rt_array;[P_vec(p_i),temp_rt_vec]];
        
    end
    
    clear MDPresult;
    
end

 %% show average neg-entropy curves as a function of sensory precision
 
figure;

p_colors = cool(length(P_vec));

for p_i = 1:length(P_vec)
    
    array2plot = H_int_array(H_int_array(:,1) == P_vec(p_i),2:end);
    array2plot = bsxfun(@minus,array2plot,array2plot(:,1));
    plot(array2plot','Color',p_colors(p_i,:));
    
%     mean_Hint = mean(H_int_array(H_int_array(:,1)==P_vec(p_i),2:end),1);
%     mean_Hint = mean_Hint - mean_Hint(1);
%     plot(mean_Hint,'Color',p_colors(p_i,:),'DisplayName',sprintf('Sensory Precision: %.2f',P_vec(p_i)));

    hold on;
end

% legend('show')



    
    
    




