%% DEM_demo_MDP_RDPsearch_LLpolicies
%  R. Conor Heins, M. Berk Mirza 2018
%
%   MDP formulation of a scene construction task (see DEM_demo_MDP_search.m) with
%   4 different scenes (hidden states), each defined as the co-occurence of 
%   two random dot patterns (RDPs) in a visual array of four quadrants, with each of the
%   two RDPs moving in orthogonal directions: ('UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP'). 
%
%    A hierarchical MDP formulation is used in order to enable independent manipulation 
%    of  uncertainty/precision at different hierarchical levels. 

% 1. The first (shallowest) hierarchical level describes the mapping from
%    a 'true' individual RDP-direction (e.g. an UP-RDP) to the observed direction. 
%    In the trivial case of no sensory uncertainty, this mapping is a bijective 
%    identity mapping, so that an RDP-direction at the higher level is directly 
%    observed as such at the lower level. The likelihood mapping at this first level,
%    from true-direction to observed direction, can however be manipulated as a model
%    of sensory uncertainty e.g. incoherent motion of RDPs.
%    In this version of RDPsearch, policies at this lower level are
%    perceptual 'decisions' about the true hidden state of the dot pattern
%    ('UP','RIGHT','DOWN','LEFT')
%
% 2. The second (deepest) hierarchical level describes the mapping from the latent scene 
%    (e.g. 'UP_RIGHT') to 12 unique manifestations of that latent scene (since there are 
%    12 different ways that an UP-RDP and a RIGHT-RDP can populate two of the four total quadrants).
%    One can manipulate the prior scene expectations (which scene does the
%    agent expect to encounter) through the variables prior_scene_belief
%    (which scene) and prior_scene_prob (percentage of trials on which agent
%    expects to see that scene).
%    The policies at this level are the categorization choices (choices
%    about the true hidden scene)

% $Id: DEM_demo_MDP_RDPsearch.m 2018-11-14 Conor Heins $
% $Id: DEM_demo_MDP_RDPsearch.m 2018-11-28 Conor Heins $
% $Id: DEM_demo_MDP_RDPsearch.m 2019-01-19 Conor Heins $

%% 1. initialize the first level: 'shallowest' MDP -- mapping RDP directions to (optionally noisy versions of) themselves in a given quadrant

% noise parameter, manipulates the coherence of dot motion via the inverse temperature parameter of a softmax transformation of
% a categorical distribution over dot directions, with a 1 at the true dot motion direction
p = 0.25; 

% temporal extent of updates at the lower level
T = 32; 

% depth of 'policies' at lower level
policy_depth = 1;

MDPshallow = initialize_MDPshallow_v6(p,T,policy_depth);
MDPshallow.chi = 1/10e20;

%% 2. initialize the second level: the'deepest' MDP -- mapping hidden scenes to RDP directions in the four quadrants

% initialize list of scene configurations and respective indices that map every 12
% configurations to one of the four scenes. 
[ all_configs,scene_idx ] = generate_scenes();

% determine agent's prior expectation over scenes, and the strength of that expectation
prior_scene_belief = 'UP_RIGHT';    % which scene does the agent expect to see
prior_scene_prob   = 0.25;          % percentage of trials that the scene is expected to occur

% choose true scene identity 
true_scene = 'UP_RIGHT';

% initialize highest level of the MDP, optionally imbuing agent with prior beliefs about
% scenes
MDPdeep = initialize_MDPdeep(prior_scene_belief,prior_scene_prob,all_configs,scene_idx,true_scene);

%% 3. Link the two levels into a single hierarchical model

% nest the shallow MDP within the deeper one
MDPdeep.MDP = MDPshallow; 
MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));

MDPfull = spm_MDP_check(MDPdeep);
clear MDPdeep MDPshallow;


%% 4. Solving

% illustrate a single trial
%==========================================================================
MDPresult = spm_MDP_VB_X(MDPfull);

%% Data analysis stuff:
% solve a sequence with a tiling of different sensory uncertainty values

% [ all_configs,scene_idx ] = generate_scenes();
% 
% prior_scene_belief = 'UP_RIGHT';    % which scene does the agent expect to see
% 
% true_scenes = {'UP_RIGHT','RIGHT_DOWN'};
% p_vec = linspace(0.05,1.25,7);
% priorprob_vec = linspace(0.25,0.8,5);
% 
% num_trials = 3;
% 
% % for debugging purposes
% true_scene_i = 1; p_i = 5; priorprob_i = 1;
% 
% for true_scene_i = 1:length(true_scenes)
%     true_scene = true_scenes{true_scene_i};
%     for p_i = 1:length(p_vec)
%         p = p_vec(p_i);
%         for priorprob_i = 1:length(priorprob_vec)
%             prior_scene_prob = priorprob_vec(priorprob_i);
%             
%             for n = 1:num_trials
% 
%                 MDPdeep = initialize_MDPdeep(prior_scene_belief,prior_scene_prob,all_configs,scene_idx,true_scene);
%                                
%                 MDPshallow = initialize_MDPshallow_witha(p);
%                 
%                 MDPdeep.MDP = MDPshallow;
%                 MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));
%                 
%                 MDPfull(1,n) = spm_MDP_check(MDPdeep);
%                 clear MDPdeep MDPshallow;
%                 
%             end
%             
%             MDPresult  = spm_MDP_VB_X(MDPfull);
%             
%             for n = 1:num_trials
%                 
%                 feedback_o = MDPresult(n).o(2,:);        
%                 results_table(true_scene_i,p_i,priorprob_i,n,strcmp(vars,'Accuracy')) =  double(any(feedback_o == 2) & ~any(feedback_o == 3));  % accuracy
%                 
%                 quadrant_o = MDPresult(n).o(1,:);
%                 quadrant_o = quadrant_o(quadrant_o > 1);
%                 num_revis = 0;
%                 uniq_quads = unique(quadrant_o);
%                 for i = 1:length(uniq_quads)
%                     quad_i_visits = length(find(quadrant_o == uniq_quads(i)));
%                     if quad_i_visits > 1
%                         num_revis = num_revis + quad_i_visits - 1;
%                     end
%                 end
%                     
%                 results_table(true_scene_i,p_i,priorprob_i,n,strcmp(vars,'Number Quadrant Revisits')) = num_revis;
%                 
%                 clear feedback_o quadrant_o num_revis uniq_quads quad_i_visits
%                 
%             end
%         end
%     end
% end
% 
           
