function MDP = DEM_demo_MDP_dot_alternative

% M. Berk Mirza
% 21-09-2018 09:25

% R. Conor Heins
% 04-10-2018 16:43

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
end

% use Gaussian distribution around center direction to parameterize sensory
% uncertainty
p = 0.5;
cert_vals = normpdf(0:3,0,p);
cert_vals = [cert_vals(1), cert_vals(2), cert_vals(2), cert_vals(3)];
cert_vals = cert_vals./sum(cert_vals);

for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)

                
        % what: A{1} -- 'n', 'U', 'R', 'D', 'L'
        %==========================================================
        if f2 == 1
            
            % at fixation location
            %----------------------------------------------------------
            A{1}(1,f1,f2) = true;
            a{1}(1,f1,f2) = power(512,10);
            
        elseif f2 > 1 && f2 < 6
            
            % saccade to cue location
            %----------------------------------------------------------
            a{1}(1,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'n'); % null
            a{1}(2,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'U'); % up
            a{1}(3,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'R'); % right
            a{1}(4,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'D'); % down
            a{1}(5,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'L'); % left
            
            if strcmp(Scenes{f1,f2 - 1},'n') == 1
                A{1}(1,f1,f2) = strcmp(Scenes{f1,f2 - 1},'n'); % null
            else
                if strcmp(Scenes{f1,f2 - 1},'U')
                    A{1}(2,f1,f2) = cert_vals(1);
                    A{1}(3,f1,f2) = cert_vals(2);
                    A{1}(5,f1,f2) = cert_vals(3);
                    A{1}(4,f1,f2) = cert_vals(4);
                end
                if strcmp(Scenes{f1,f2-1},'R')
                    A{1}(3,f1,f2) = cert_vals(1);
                    A{1}(2,f1,f2) = cert_vals(2);
                    A{1}(4,f1,f2) = cert_vals(3);
                    A{1}(5,f1,f2) = cert_vals(4);
                end
                if strcmp(Scenes{f1,f2-1},'D')
                    A{1}(4,f1,f2) = cert_vals(1);
                    A{1}(3,f1,f2) = cert_vals(2);
                    A{1}(5,f1,f2) = cert_vals(3);
                    A{1}(2,f1,f2) = cert_vals(4);
                end
                if strcmp(Scenes{f1,f2-1},'L')
                    A{1}(5,f1,f2) = cert_vals(1);
                    A{1}(4,f1,f2) = cert_vals(2);
                    A{1}(2,f1,f2) = cert_vals(3);
                    A{1}(3,f1,f2) = cert_vals(4);
                end
            end
                       
        elseif f2 > 5
            
            % saccade choice location
            %------------------------------------------------------
            A{1}(6,f1,f2) = ceil(f1/num_configs) + 5 == f2; % right
            A{1}(7,f1,f2) = ceil(f1/num_configs) + 5 ~= f2; % wrong
            
            a{1}(6,f1,f2) = power(512,10)*(ceil(f1/num_configs) + 5 == f2); % right
            a{1}(7,f1,f2) = power(512,10)*(ceil(f1/num_configs) + 5 ~= f2); % wrong
            
        end
        
        % where: A{2} {'start','1',...,'4','UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP'}
        %----------------------------------------------------------
        A{2}(f2,f1,f2) = 1; % by making this always 1, you make it such that hidden states of motor plants (eye positions) are directly observed
        % i.e. there's no uncertainty about where you're looking
        
        a{2}(f2,f1,f2) = power(512,10);

    end
end
for g = 1:Ng
    A{g} = double(A{g});
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
mdp.T = T;                      % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.a = a;                      % likelihood mapping in the generative model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = d;                      % prior over initial states
mdp.s = [1 1]';                 % initial states
mdp.o = [1 1]';                 % initial outcome

mdp.Aname = {'what','where'};
mdp.Bname = {'category_configuration','where'};
 
mdp = spm_MDP_check(mdp);

%% illustrate a sequence of trials
%==========================================================================
 
% true initial states –
%==========================================================================

tic
MDP = spm_MDP_VB_X(mdp);
toc

%% solve and evaluate performance for sequence of trials

N = 500;
s(1,:) = ceil(rand(1,N)*48);

for i = 1:N
    MDP(i)   = mdp;      % create structure array
    MDP(i).s(1) = s(i);   % context
end

MDP = spm_MDP_VB_X(MDP);


% illustrate a single trial
%==========================================================================
MDP   = spm_MDP_VB_X(mdp);
 
% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP);
 
% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],1);

return
 