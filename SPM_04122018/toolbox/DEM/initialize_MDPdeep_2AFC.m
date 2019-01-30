function [ MDPdeep ] = initialize_MDPdeep_2AFC(true_scene)
%INITIALIZE_MDPDEEP Initialize deep level of MDP, and equip with a true
%latent state (true_scene)
% INPUTS:  true_scene:        true scene identity ('LEFT' or 'RIGHT')
% OUTPUTS: MDPdeep:           initialized higher/deeper level of hierarchical MDP


switch true_scene
    case 'LEFT'
        scene_idx = 1;
    case 'RIGHT'
        scene_idx = 2;
end

D{1} = ones(2,1);                         % what:     True hidden state of the world: options are {'LEFT','RIGHT'}
D{2} = [1 0 0 0]';                    % where:    {'fixation-start','Target Area','LEFT choice','RIGHT choice'}

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end

% outcome dimensions 
% first visual dimension: {'Null', 'LEFT','RIGHT'}; 3
% second trial-feedback dimension: {'No Choice Yet','Correct','Incorrect'} 2
% third proprioceptive dimension: {'Fixation-start','Target Area','LEFT
% choice', 'RIGHT choice'} 4

No    = [3 3 4]; 
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
            A{2}(1,f1,f2)   = 1; 
            
        elseif f2 > 1 && f2 < 3
            % saccade to cue location
            %----------------------------------------------------------

            if f1 == 1 % in the case that the true scene is 'LEFT'
                A{1}(2,f1,f2)   = 1;
            elseif f1 == 2 % in the case that the true scene is 'RIGHT'
                A{1}(3,f1,f2) = 1;
            end
            
            A{2}(1,f1,f2) = true; % if you haven't looked at the choice area yet, you're undecided
            
        elseif f2 >= 3
            
            % saccade choice location
            %------------------------------------------------------
            
            if f1 == 1 % if true scene is LEFT
                if f2 == 3
                    A{2}(2,f1,f2) = 1; % if you're looking at 'LEFT CHOICE' (i.e. f2 == 3), then you're correct
                else
                    A{2}(3,f1,f2) = 1; % if you're looking anywhere else than 'LEFT CHOICE,' you're wrong (hence *FORCED* part of 2-AFC)
                end
            elseif f1 == 2 % if true scene is RIGHT
                if f2 == 4
                    A{2}(2,f1,f2) = 1; % if you're looking at 'RIGHT CHOICE' (i.e. f2 == 4), then you're correct
                else
                    A{2}(3,f1,f2) = 1; % if you're looking anywhere else than 'RIGHT CHOICE,' you're wrong (hence *FORCED* part of 2-AFC)
                end
            end

            A{1}(1,f1,f2)   = 1;       % when you're looking at either of the choice areas, your visual feedback is null (same as looking at fixation center)
            
        end
        
        % where: A{2} {'Fixation-start','Target Area','LEFT choice', 'RIGHT choice'}
        %----------------------------------------------------------
        A{3}(f2,f1,f2) = 1; % you're proprioceptive feedback is always identity
          
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
T         = 3; % they have to choose after only one fixation on the target area
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{3}      = zeros(No(3),T);
C{2}(1,:) =  0;                 % the agent doesn't expect to be undecided (can make this 0 or negative, try different values)
C{2}(2,:) =  2;                 % the agent expects to be right
C{2}(3,:) = -4;                 % the agent would be really surprised to be wrong


% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
MDPdeep.T = T;                      % number of moves
MDPdeep.U = U;                      % allowable policies
MDPdeep.A = A;                      % observation model
MDPdeep.B = B;                      % transition probabilities
MDPdeep.C = C;                      % preferred outcomes
MDPdeep.D = D;                      % prior over initial states
MDPdeep.s = [scene_idx 1]';        % initial states

MDPdeep.Aname = {'what','feedback','where'};
MDPdeep.Bname = {'category_configuration','where'};
 
MDPdeep = spm_MDP_check(MDPdeep);

end

