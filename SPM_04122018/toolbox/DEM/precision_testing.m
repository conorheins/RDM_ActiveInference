%% 1. initialize the first level: 'shallowest' MDP -- mapping RDP directions to (optionally noisy versions of) themselves in a given quadrant

p = 0.5; % noise parameter, standard deviation of Gaussian distribution centered around RDP motion direction

MDPshallow = initialize_MDPshallow_witha(p);

%% 2. initialize the second level: the'deepest' MDP -- mapping hidden scenes to RDP directions in the four quadrants

% initialize list of scene configurations and respective indices that map every 12
% configurations to one of the four scenes. 
[ all_configs,scene_idx ] = generate_scenes();

% determine agent's prior expectation over scenes, and the strength of that expectation
prior_scene_belief = 'RIGHT_DOWN';    % which scene does the agent expect to see
prior_scene_prob   = 0.6;           % percentage of trials that the scene is expected to occur

% now repeat but changing prior precision
%--------------------------------------------------------------------------
b     = (0:8)/8;
b     = 1./(2*b + 1/8);

% choose true scene identity 
true_scene = 'UP_RIGHT';

all_precisions = cell(1,length(b));

for j = 1:length(b)
    
    for i = 1:10
        
        MDPdeep = initialize_MDPdeep(prior_scene_belief,prior_scene_prob,all_configs,scene_idx,true_scene);
        MDPdeep.beta = b(j);
        MDPdeep.MDP = MDPshallow; 
        MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));
        
        MDPfull(1,i) = spm_MDP_check(MDPdeep);
        clear MDPdeep;
        
    end
    
    all_precisions{1,j}  = spm_MDP_VB_X(MDPfull); clear MDPfull;
    
end

%% look at some of these trials, mayne!

for j = 2:length(b)
    
    this_precis = all_precisions{1,j};
    
    for i = 1:size(this_precis,2)
        fprintf('Trial %d of Precision value: %.2f\n',i,1/b(j));
        fprintf('Prior Belief: %s, Strength: %.2f, Sensory Uncert: %.2f\n',...
            prior_scene_belief,prior_scene_prob,p);
        MDP_beliefs_video(this_precis,i,all_configs);
        pause;
        close gcf;
    end
    
end

    