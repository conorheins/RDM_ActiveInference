function [mdp,Scenes] = initialize_searchMDP_4(p,s,prior_scene,prior_scene_prob)

%% Initialize prior beliefs about initial states (in terms of counts): D and d
%--------------------------------------------------------------------------

num_cats = 4; num_configs = 12;
d{1} = ones(num_cats*num_configs,1);            % what:     Categories {'UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP'} and Configurations (there are twelve unique configurations)
d{2} = [1 0 0 0 0 0 0 0 0]';                    % where:    {'start','1',...,'4','UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP'}

U_R = unique(perms([2 3 1 1]),'rows');
R_D = unique(perms([3 4 1 1]),'rows');
D_L = unique(perms([4 5 1 1]),'rows');
L_U = unique(perms([5 2 1 1]),'rows');

Scene = zeros(num_cats*num_configs,size(U_R,2));
unique_cat_configs = [1:num_configs:(num_configs*num_cats),num_configs*num_cats];
Scene(unique_cat_configs(1):(unique_cat_configs(2)-1),:) = U_R;
Scene(unique_cat_configs(2):(unique_cat_configs(3)-1),:) = R_D;
Scene(unique_cat_configs(3):(unique_cat_configs(4)-1),:) = D_L;
Scene(unique_cat_configs(4):unique_cat_configs(5),:) = L_U;

d_val = ((prior_scene_prob * num_cats - prior_scene_prob)/(1-prior_scene_prob));
switch prior_scene
    case 'UP_RIGHT'
        d{1}(unique_cat_configs(1):(unique_cat_configs(2)-1)) = d_val;
    case 'RIGHT_DOWN'
        d{1}(unique_cat_configs(2):(unique_cat_configs(3)-1)) = d_val;
    case 'DOWN_LEFT'
        d{1}(unique_cat_configs(3):(unique_cat_configs(4)-1)) = d_val;
    case 'LEFT_UP'
        d{1}(unique_cat_configs(4):unique_cat_configs(5)) = d_val;
end
        

Scenes = num2cell(Scene);
Scenes(find([Scenes{:}] == 1)) = {'n'}; 
Scenes(find([Scenes{:}] == 2)) = {'U'};
Scenes(find([Scenes{:}] == 3)) = {'R'};
Scenes(find([Scenes{:}] == 4)) = {'D'};
Scenes(find([Scenes{:}] == 5)) = {'L'};

%% Generate likelihood matrix A: probabilistic mapping from hidden states to outcomes
%--------------------------------------------------------------------------
Nf    = numel(d);
for f = 1:Nf
    Ns(f) = numel(d{f});
end
No    = [7 9];
Ng    = numel(No);
for g = 1:Ng
    A{g} = zeros([No(g),Ns]);
    GP{g} = A{g};
end

% use Gaussian distribution around center direction to parameterize sensory
% uncertainty
if p == 0
    cert_vals = [1 0 0 0 ];
else
    cert_vals = normpdf(0:3,0,p);
    cert_vals = [cert_vals(1), cert_vals(2), cert_vals(2), cert_vals(3)];
    cert_vals = cert_vals./sum(cert_vals);
end

for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)

                
        % what: A{1} -- 'n', 'U', 'R', 'D', 'L'
        %==========================================================
        if f2 == 1
            
            % at fixation location
            %----------------------------------------------------------
            A{1}(1,f1,f2) = true;
            
        elseif f2 > 1 && f2 < 6
            
            % saccade to cue location
            %----------------------------------------------------------
            A{1}(1,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'n'); % null
            A{1}(2,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'U'); % up
            A{1}(3,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'R'); % right
            A{1}(4,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'D'); % down
            A{1}(5,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'L'); % left
            
            if strcmp(Scenes{f1,f2 - 1},'n') == 1
                GP{1}(1,f1,f2) = strcmp(Scenes{f1,f2 - 1},'n'); % null
            else
                if strcmp(Scenes{f1,f2 - 1},'U')
                    GP{1}(2,f1,f2) = cert_vals(1);
                    GP{1}(3,f1,f2) = cert_vals(2);
                    GP{1}(5,f1,f2) = cert_vals(3);
                    GP{1}(4,f1,f2) = cert_vals(4);
                end
                if strcmp(Scenes{f1,f2-1},'R')
                    GP{1}(3,f1,f2) = cert_vals(1);
                    GP{1}(2,f1,f2) = cert_vals(2);
                    GP{1}(4,f1,f2) = cert_vals(3);
                    GP{1}(5,f1,f2) = cert_vals(4);
                end
                if strcmp(Scenes{f1,f2-1},'D')
                    GP{1}(4,f1,f2) = cert_vals(1);
                    GP{1}(3,f1,f2) = cert_vals(2);
                    GP{1}(5,f1,f2) = cert_vals(3);
                    GP{1}(2,f1,f2) = cert_vals(4);
                end
                if strcmp(Scenes{f1,f2-1},'L')
                    GP{1}(5,f1,f2) = cert_vals(1);
                    GP{1}(4,f1,f2) = cert_vals(2);
                    GP{1}(2,f1,f2) = cert_vals(3);
                    GP{1}(3,f1,f2) = cert_vals(4);
                end
            end
                       
        elseif f2 > 5
            
            % saccade choice location
            %------------------------------------------------------
            GP{1}(6,f1,f2) = ceil(f1/num_configs) + 5 == f2; % right
            GP{1}(7,f1,f2) = ceil(f1/num_configs) + 5 ~= f2; % wrong
            
            A{1}(6,f1,f2) = (ceil(f1/num_configs) + 5 == f2); % right
            A{1}(7,f1,f2) = (ceil(f1/num_configs) + 5 ~= f2); % wrong
            
        end
        
        % where: A{2} {'start','1',...,'4','UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP'}
        %----------------------------------------------------------
        GP{2}(f2,f1,f2) = 1; 
        A{2}(f2,f1,f2) = true;

    end
end
for g = 1:Ng
    A{g} = double(A{g});
    GP{g} = double(GP{g});
end
 
%% controlled transitions: B{f} for each factor
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
 
%% allowable policies (here, specified as the next action) U
%--------------------------------------------------------------------------
Np        = size(B{2},3);
U         = ones(1,Np,Nf);
U(:,:,2)  = 1:Np;
 
% priors: (utility) C
%--------------------------------------------------------------------------
T         = 8;
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{1}(6,:) =  2;                 % the agent expects to be right
C{1}(7,:) = -4;                 % and not wrong

% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.T  = T;                      % number of moves
mdp.U  = U;                      % allowable policies
mdp.A  = A;                      % observation model
mdp.GP = GP;                     % generative process
mdp.B  = B;                      % transition probabilities
mdp.C  = C;                      % preferred outcomes
mdp.D  = d;                      % prior over initial states
mdp.s  = [s 1]';                 % initial states
o      = [1 1]';                 % initial outcome, in terms of indices

mdp.Aname = {'what','where'};
mdp.Bname = {'category_configuration','where'};
 
mdp = spm_MDP_check(mdp);

end