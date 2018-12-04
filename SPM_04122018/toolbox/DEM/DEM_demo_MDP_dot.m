function MDP = DEM_demo_MDP_dot

%% prior beliefs about initial states (in terms of counts_: D and d
%--------------------------------------------------------------------------
d{1} = [1 1 1 1]';            % what:     {'A','B','C','D'}
d{2} = ones(12,1);            % Scene:    Configurations (there are twelve unique configurations)
d{3} = [1 0 0 0 0 0 0 0 0]';  % where:    {'start','1',...,'4','A','B','C','D'}

sA = num2cell(unique(perms([1 2 3 3]),'rows'));
sA(find([sA{:}] == 1)) = {'u'};
sA(find([sA{:}] == 2)) = {'r'};
sA(find([sA{:}] == 3)) = {'n'};

sB = num2cell(unique(perms([1 2 3 3]),'rows'));
sB(find([sB{:}] == 1)) = {'r'};
sB(find([sB{:}] == 2)) = {'d'};
sB(find([sB{:}] == 3)) = {'n'};

sC = num2cell(unique(perms([1 2 3 3]),'rows'));
sC(find([sC{:}] == 1)) = {'d'};
sC(find([sC{:}] == 2)) = {'l'};
sC(find([sC{:}] == 3)) = {'n'};

sD = num2cell(unique(perms([1 2 3 3]),'rows'));
sD(find([sD{:}] == 1)) = {'l'};
sD(find([sD{:}] == 2)) = {'u'};
sD(find([sD{:}] == 3)) = {'n'};

Scenes{1} = sA;
Scenes{2} = sB;
Scenes{3} = sC;
Scenes{4} = sD;

%% probabilistic mapping from hidden states to outcomes: A
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
        for f3 = 1:Ns(3)
                
            % what: A{1}
            %==========================================================
            if f3 == 1

                % at fixation location
                %----------------------------------------------------------
                A{1}(1,f1,f2,f3) = true;

            elseif f3 > 1 && f3 < 6

                % saccade to cue location
                %----------------------------------------------------------
                A{1}(1,f1,f2,f3) = strcmp(Scenes{f1}{f2,f3 - 1},'n'); % null
                A{1}(2,f1,f2,f3) = strcmp(Scenes{f1}{f2,f3 - 1},'u'); % up
                A{1}(3,f1,f2,f3) = strcmp(Scenes{f1}{f2,f3 - 1},'r'); % right
                A{1}(4,f1,f2,f3) = strcmp(Scenes{f1}{f2,f3 - 1},'d'); % down
                A{1}(5,f1,f2,f3) = strcmp(Scenes{f1}{f2,f3 - 1},'l'); % left

            elseif f3 > 5

                % saccade choice location
                %------------------------------------------------------
                A{1}(6,f1,f2,f3) = (f3 - 5) == f1; % right
                A{1}(7,f1,f2,f3) = (f3 - 5) ~= f1; % wrong

            end

            % where: A{2} {'start','1',...,'4','A','B','C','D'}
            %----------------------------------------------------------
            A{2}(f3,f1,f2,f3) = 1;
                
        end
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

%% controllable fixation points: move to the k-th location
%--------------------------------------------------------------------------
for k = 1:Ns(3)
    B{3}(:,:,k) = 0;
    B{3}(k,:,k) = 1;
end
 
%% allowable policies (here, specified as the next action) U
%--------------------------------------------------------------------------
Np        = size(B{3},3);
U         = ones(1,Np,Nf);
U(:,:,3)  = 1:Np;
 
%% priors: (utility) C
%--------------------------------------------------------------------------
T         = 6;
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{1}(6,:) =  2;                 % the agent expects to be right
C{1}(7,:) = -4;                 % and not wrong
C{2}(1:5,5:end) = -4;           % make tardy sampling costly


%% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.T = T;                      % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = d;                      % prior over initial states
mdp.s = [1 1 1]';               % initial state (if state is e.g. [1, 1, 1], this means: Scene A, Configuration 1, Eyes start at Fixation)
mdp.o = [1 1]';                 % initial outcome

mdp.Aname = {'what','where'};
mdp.Bname = {'what','scene','where'};
 
mdp = spm_MDP_check(mdp);
%% illustrate a single trial
%==========================================================================
MDP   = spm_MDP_VB_X_old(mdp);

% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP);

view_trial_trajectory(MDP,1,2);
 
% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],1);

% illustrate a sequence of trials
%==========================================================================
 
% true initial states –
%==========================================================================
N = 10000;
s(1,:) = ceil(rand(1,N)*4);
s(2,:) = ceil(rand(1,N)*12);

for i = 1:N
    MDP(i)   = mdp;      % create structure array
    MDP(i).s(1:2,1) = s(:,i);   % context
end

% Solve - an example sequence
%==========================================================================
MDP  = spm_MDP_VB_X_old(MDP);

save('MDPresults_09192018.mat','MDP');

scene_names = {'SceneA', 'SceneB', 'SceneC', 'SceneD'};
for i = 1:N
    view_trial_trajectory(MDP,i,2);
    pause;
    close gcf;
end
% 
% scene_names = {'SceneA', 'SceneB', 'SceneC', 'SceneD'};
% for i = 1:N
%     view_trial_trajectory(MDP,i,2);
%     states = MDP(i).s(:,1);
%     saveas(gcf,fullfile('trial_trajectories_09182018',sprintf('%s_config%d_trial%d.png',scene_names{states(1)},states(2),i)));
%     close gcf;
% end

% solve and evaluate performance
%----------------------------------------------------------------------
for i = 1:N
    o      = MDP(i).o(1,:);                       % outcomes
    P(1,i) = double(any(o == 6) & ~any(o == 7));  % accuracy
    P(2,i) = find([(o > 5), 1],1) - 1;            % number of saccades
    P(3,i) = mean(MDP(i).rt);                     % reaction time
end

result_table = zeros(1000,5);
for i = 1:N
    result_table(i,1) = MDP(i).s(1,1);
    result_table(i,2) = MDP(i).s(2,1);
    o      = MDP(i).o(1,:);                       % outcomes
    result_table(i,3) = double(any(o == 6) & ~any(o == 7));  % accuracy
    result_table(i,4) = find([(o > 5), 1],1) - 1;            % number of saccades
    result_table(i,5) = mean(MDP(i).rt);                     % reaction time
end


% solve and evaluate performance
%----------------------------------------------------------------------
averages = mean(P,2);

return
 