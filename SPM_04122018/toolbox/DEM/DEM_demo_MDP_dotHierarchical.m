function MDP = DEM_demo_MDP_dot

% Conor Heins, M. Berk Mirza

% first level
% Observe the dot motion. There is one hidden state dimension that 
% corresponds to dot direction. This directly maps onto the dot direction 
% outcome.

p = 0;
% use Gaussian distribution around center direction to parameterize sensory
% uncertainty
if p == 0
    cert_vals = [1 0 0 0];
else
    cert_vals = normpdf(0:3,0,p);
    cert_vals = [cert_vals(1), cert_vals(2), cert_vals(2), cert_vals(3)];
    cert_vals = cert_vals./sum(cert_vals);
end

D{1} = [ 1 1 1 1 1 ]'; % Direction: {'Null','Up','Right','Down','Left'}

Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end

for f1 = 1:Ns(1)
    if f1 == 1
        A{1}(f1,f1) = 1;
    elseif f1 == 2 % U
        A{1}(2,f1) = cert_vals(1); % U
        A{1}(4,f1) = cert_vals(4); % D
        A{1}(3,f1) = cert_vals(2); % R
        A{1}(5,f1) = cert_vals(3); % L
    elseif f1 == 3 % R
        A{1}(3,f1) = cert_vals(1); % R
        A{1}(5,f1) = cert_vals(4); % L
        A{1}(2,f1) = cert_vals(2); % U
        A{1}(4,f1) = cert_vals(3); % D
    elseif f1 == 4 % D
        A{1}(4,f1) = cert_vals(1); % D
        A{1}(2,f1) = cert_vals(4); % U
        A{1}(3,f1) = cert_vals(2); % R
        A{1}(5,f1) = cert_vals(3); % L
    elseif f1 == 5 % L
        A{1}(5,f1) = cert_vals(1); % L
        A{1}(3,f1) = cert_vals(4); % R
        A{1}(2,f1) = cert_vals(2); % U
        A{1}(4,f1) = cert_vals(3); % D
    end 
end

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% MDP Structure
%--------------------------------------------------------------------------
mdp.T = 1;                      % number of updates
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.D = D;                      % prior over initial states

mdp.Aname = {'Direction'};
mdp.Bname = {'Direction'};
clear A B D
 
MDP = spm_MDP_check(mdp);
clear mdp

% Second level
% Infer the scene category based on observations
% There are two hidden state dimensions. The first is the scene category 
% and scene configuration. The second corresponds to the locations in the
% scene. These hidden state dimensions determine the dot directions U, R, D
% and L. This translates to the initial states in the subordinate level.
% There are three outcome modalities. First is the dot direction, second is
% a feedback about whether the scene is categorised correctly or not. This
% outcome modality has three outcomes (null, right and wrong). This stays
% as null until the agent makes a decision about the scene category. The 
% third is the where modality.

% prior beliefs about initial states (in terms of counts): D and d
%--------------------------------------------------------------------------
D{1} = ones(48,1);            % what:     {'A','B','C','D'} and Scene Configurations (there are twelve unique configurations)
D{2} = [1 0 0 0 0 0 0 0 0]';  % where:    {'start','1',...,'4','A','B','C','D'}

sA = unique(perms([2 3 1 1]),'rows');
sB = unique(perms([3 4 1 1]),'rows');
sC = unique(perms([4 5 1 1]),'rows');
sD = unique(perms([5 2 1 1]),'rows');

Scene = zeros(size(D{1},1),size(sA,2));
Scene(1:12,:)  = sA; Scene(13:24,:) = sB; Scene(25:36,:) = sC; Scene(37:48,:) = sD;
Scenes = num2cell(Scene);
Scenes(find([Scenes{:}] == 1)) = {'n'}; 
Scenes(find([Scenes{:}] == 2)) = {'u'};
Scenes(find([Scenes{:}] == 3)) = {'r'};
Scenes(find([Scenes{:}] == 4)) = {'d'};
Scenes(find([Scenes{:}] == 5)) = {'l'};

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
            
            A{1}(1,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'n'); % null
            A{1}(2,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'u'); % up
            A{1}(3,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'r'); % right
            A{1}(4,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'d'); % down
            A{1}(5,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'l'); % left
            A{2}(1,f1,f2)   = 1; % null
            
        elseif f2 > 5
            
            % saccade choice location
            %------------------------------------------------------
            A{2}(2,f1,f2) = ceil(f1/12) + 5 == f2; % right
            A{2}(3,f1,f2) = ceil(f1/12) + 5 ~= f2; % wrong
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

% MDP Structure
%--------------------------------------------------------------------------
mdp.MDP  = MDP;
mdp.link = sparse(1,1,1,numel(MDP.D),Ng);

% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.T = T;                      % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states
mdp.s = [1 1]';                 % initial states

mdp.Aname = {'what','feedback','where'};
mdp.Bname = {'category-configuration','where'};
 
mdp = spm_MDP_check(mdp);
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
 