%% spm_MDP_trust_network
%  Brennan Klein, Conor Heins 2018
%
% This routine uses the Markov decision process formulation of active
% inference (with variational Bayes) to model a trust game. In trust
% games, one plays an opponent who can either cooperate or defect. The
% payoff contingencies depend upon the joint choices of you and your
% opponent, which in turn depend upon your inferences about the nature of
% the opponent (pro-social or non-social). This example illustrates
% a network of agents playing these trust games and extends the 
% original formulation (see 'spm_MDP_trust.m') to include a deep
% hierarchical model, wherein the outcomes of lower level games
% affect the agent's latent beliefs about the presence of
% pro- vs. anti-social in their network environment.
%
%    In the following, we detail the nature of states/outcomes at the two
%    hierarchical levels of the agents' generative models
%
% 1. The first (shallowest) hierarchical level describes the mapping from
%    9 hidden states to 5 'reward' outcomes (cooperate-cooperate, cooperate-defect, 
%    defect-cooperate, and defect-defect -- describing the combination of
%    the agent and its opponent's choices). Each outcome is associated with
%    a prior expectation (or preference), that is derived in the classical
%    way from the game's payoff matrix. The first hidden state is a starting state 
%    and the subsequent eight states model the four combinations of cooperation and defection
%    (between you and your opponent), each combined with a prior belief that the opponent 
%    is either pro-social or non-social. Initially, these prior beliefs are uninformative 
%    but are subsequently informed through experience and the states of the 'hyper-prior'
%    distribution over the general sociality of agents, which is described at the deeper
%    (second) level of the hierarchical model.
% 2. The second (deepest) hierarchical level describes the mapping from
%    latent beliefs about the expected 'pro- vs. anti-socialness' of agents encountered
%    in each agent's local network environment, and the posterior expectations about the hidden states 
%    of the world, which then become  the hidden states for the lower level
%    of the hierarchical model


%% 1. initialize the first level: 'shallowest' MDP -- mapping 9 hidden states to individual game outcomes

k = 2; % number of opponents (number of 'edges' or 'neighbors' in the network)

U{1} = [26 10;                      % self payoff (utility)
        21 18]/8;                   
        
U{2} = [26 42;                      % other payoff (utility)
         7 10]/8;                   

MDPshallow = initialize_trustMDP_shallow(k,U); 

%% 2. initialize the second level: the'deepest' MDP -- mapping hidden scenes to RDP directions in the four quadrants

% initialize highest level of the MDP, optionally imbuing agent with prior beliefs about
% scenes
MDPdeep = initialize_MDPdeep(prior_scene_belief,prior_scene_prob,all_configs,scene_idx,true_scene);

%% 3. Link the two levels into a single hierarchical model

% nest the shallow MDP within the deeper one
MDPdeep.MDP = MDPshallow; 
MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));

MDPfull = spm_MDP_check(MDPdeep);

%% 4. Solving

% illustrate a single trial
%==========================================================================
MDPresult  = spm_MDP_game(MDPfull);