function MDP = DEM_demo_MDP_dot

% M. Berk Mirza
% 21-09-2018 09:25

% this term controls how coherently the dots will move (p=0 means
% completely incoherently, p = 16 means completely coherently)
p = 16;

% prior beliefs about initial states (in terms of counts): D and d
%--------------------------------------------------------------------------
d{1} = ones(48,1);            % what:     {'A','B','C','D'} and Scene Configurations (there are twelve unique configurations)
d{2} = [1 0 0 0 0 0 0 0 0]';  % where:    {'start','1',...,'4','A','B','C','D'}

sA = unique(perms([2 3 1 1]),'rows');
sB = unique(perms([3 4 1 1]),'rows');
sC = unique(perms([4 5 1 1]),'rows');
sD = unique(perms([5 2 1 1]),'rows');

Scene = zeros(size(d{1},1),size(sA,2));
Scene(1:12,:)  = sA; Scene(13:24,:) = sB; Scene(25:36,:) = sC; Scene(37:48,:) = sD;
Scenes = num2cell(Scene);
Scenes(find([Scenes{:}] == 1)) = {'n'}; 
Scenes(find([Scenes{:}] == 2)) = {'u'};
Scenes(find([Scenes{:}] == 3)) = {'r'};
Scenes(find([Scenes{:}] == 4)) = {'d'};
Scenes(find([Scenes{:}] == 5)) = {'l'};

% probabilistic mapping from hidden states to outcomes: A
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
for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)
                
        % what: A{1}
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
            a{1}(2,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'u'); % up
            a{1}(3,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'r'); % right
            a{1}(4,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'d'); % down
            a{1}(5,f1,f2)   = power(512,10)*strcmp(Scenes{f1,f2 - 1},'l'); % left
            
            if strcmp(Scenes{f1,f2 - 1},'n') == 1
                A{1}(1,f1,f2) = strcmp(Scenes{f1,f2 - 1},'n'); % null
            else
                A{1}(2,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'u'); % up
                A{1}(3,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'r'); % right
                A{1}(4,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'d'); % down
                A{1}(5,f1,f2)   = strcmp(Scenes{f1,f2 - 1},'l'); % left
                A{1}(2:5,f1,f2) = spm_softmax(A{1}(2:5,f1,f2)*p); 
            end
            
        elseif f2 > 5
            
            % saccade choice location
            %------------------------------------------------------
            A{1}(6,f1,f2) = ceil(f1/12) + 5 == f2; % right
            A{1}(7,f1,f2) = ceil(f1/12) + 5 ~= f2; % wrong
            
            a{1}(6,f1,f2) = power(512,10)*(ceil(f1/12) + 5 == f2); % right
            a{1}(7,f1,f2) = power(512,10)*(ceil(f1/12) + 5 ~= f2); % wrong
            
        end
        
        % where: A{2} {'start','1',...,'4','A','B','C','D'}
        %----------------------------------------------------------
        A{2}(f2,f1,f2) = 1;
        a{2}(f2,f1,f2) = power(512,10);
                
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
 