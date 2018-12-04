function [ MDPdeep ] = initialize_MDPdeep( prior_scene_belief,prior_scene_prob,all_configs,scene_idx,true_scene)
%INITIALIZE_MDPDEEP Initialize deep level of MDP, optionally equipping the agent
% with prior expectations about scene occurrence probability (through the
% parameters prior_scene_belief and prior_scene_prob).
% INPUTS: prior_scene_belief: which scene does the agent believe
%                             predominates trials
%                             ('UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP')
%         prior_scene_prob:   how likely  does the agent believe
%                             'prior_scene_belief' is (percentage, between
%                             0 and 1)
%         all_configs:        cell array of all unique configurations, in terms of what RDP ('n','U','R','D','L') occupies 
%                             each of the four quadrants for a given configuration (48 x 4 cell array)
%         scene_idx:          cell array of indices that map 12 unique
%                             configurations in all_configs to one hidden scene
%         true_scene:         true scene identity ('UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP')
%
% OUTPUTS: MDPdeep:           initialized higher/deeper level of hierarchical MDP


num_scenes = length(scene_idx); 
num_configs = size(all_configs,1)/num_scenes; 

if ~exist('prior_scene_belief','var') || isempty(prior_scene_belief)
    prior_scene_belief = 'UP_RIGHT';
end

if ~exist('prior_scene_prob','var') || isempty(prior_scene_prob)
    prior_scene_prob = 1/num_scenes; % in absence of prior probability, defaults to 'flat' prior over scene identity 
end

if ~exist('all_configs','var') || ~exist('scene_idx','var')
    error('Please use generate_scenes() to create variables all_configs and scene_idx')
end

if ~exist('true_scene','var') || isempty(true_scene)
    true_scene = ceil(rand(1)*num_scenes*num_configs);
elseif ischar(true_scene)
    switch true_scene
        case 'UP_RIGHT'
            true_scene = randi([scene_idx{1}(1) scene_idx{1}(end)],1);
        case 'RIGHT_DOWN'
            true_scene = randi([scene_idx{2}(1) scene_idx{2}(end)],1);
        case 'DOWN_LEFT'            
            true_scene = randi([scene_idx{3}(1) scene_idx{3}(end)],1);
        case 'LEFT_UP'
            true_scene = randi([scene_idx{4}(1) scene_idx{4}(end)],1);
    end
end


D{1} = ones(num_scenes*num_configs,1);             % what:     Categories {'UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP'} and Configurations (there are twelve unique configurations)
D{2} = [1 0 0 0 0 0 0 0 0]';                    % where:    {'start','1',...,'4','UP_RIGHT choice','RIGHT_DOWN choice','DOWN_LEFT choice','LEFT_UP choice'}

d_belief = ((prior_scene_prob * num_scenes - prior_scene_prob)/(1-prior_scene_prob));
switch prior_scene_belief
    case 'UP_RIGHT'
        D{1}(scene_idx{1}) = d_belief;
    case 'RIGHT_DOWN'
        D{1}(scene_idx{2}) = d_belief;
    case 'DOWN_LEFT'
        D{1}(scene_idx{3}) = d_belief;
    case 'LEFT_UP'
        D{1}(scene_idx{4}) = d_belief;
end

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end
No    = [5 3 9];
Ng    = numel(No);
for g = 1:Ng
    A{g} = zeros([No(g),Ns]);
end
for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)
                
        % what: A{1}
        %==========================================================
        if f2 == 1
            
            % at fixation location
            %----------------------------------------------------------
            A{1}(1,f1,f2) = true;
            A{2}(1,f1,f2)   = 1; % null
            
        elseif f2 > 1 && f2 < 6
            
            % saccade to cue location
            %----------------------------------------------------------
            
            A{1}(1,f1,f2)   = strcmp(all_configs{f1,f2 - 1},'n'); % null
            A{1}(2,f1,f2)   = strcmp(all_configs{f1,f2 - 1},'U'); % UP
            A{1}(3,f1,f2)   = strcmp(all_configs{f1,f2 - 1},'R'); % RIGHT
            A{1}(4,f1,f2)   = strcmp(all_configs{f1,f2 - 1},'D'); % DOWN
            A{1}(5,f1,f2)   = strcmp(all_configs{f1,f2 - 1},'L'); % LEFT
            A{2}(1,f1,f2)   = 1;                             % null
            
        elseif f2 > 5
            
            % saccade choice location
            %------------------------------------------------------
            A{2}(2,f1,f2) = ceil(f1/num_configs) + 5 == f2; % right
            A{2}(3,f1,f2) = ceil(f1/num_configs) + 5 ~= f2; % wrong
            A{1}(1,f1,f2)   = 1; % null
            
        end
        
        % where: A{2} {'start','1',...,'4','A','B','C','D'}
        %----------------------------------------------------------
        A{3}(f2,f1,f2) = 1;
          
    end
end
for g = 1:Ng
    A{g} = double(A{g});
end

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% controllable fixation points: move to the k-th location
%--------------------------------------------------------------------------
for k = 1:Ns(2)
    B{2}(:,:,k) = 0;
    B{2}(k,:,k) = 1;
end
 
% allowable policies (here, specified as the next action) U
%--------------------------------------------------------------------------
Np        = size(B{2},3);
U         = ones(1,Np,Nf);
U(:,:,2)  = 1:Np;
 
% priors: (utility) C
%--------------------------------------------------------------------------
T         = 8;
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{3}      = zeros(No(3),T);
C{2}(2,:) =  2;                 % the agent expects to be right
C{2}(3,:) = -4;                 % and not wrong


% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
MDPdeep.T = T;                      % number of moves
MDPdeep.U = U;                      % allowable policies
MDPdeep.A = A;                      % observation model
MDPdeep.B = B;                      % transition probabilities
MDPdeep.C = C;                      % preferred outcomes
MDPdeep.D = D;                      % prior over initial states
MDPdeep.s = [true_scene 1]';                 % initial states

MDPdeep.Aname = {'what','feedback','where'};
MDPdeep.Bname = {'category_configuration','where'};
 
MDPdeep = spm_MDP_check(MDPdeep);

end

