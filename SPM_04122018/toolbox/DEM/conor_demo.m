%% DEM_demo_conortest

% prior beliefs over states 
% three states of the world (Goo, Gah, and Gi)
% four control states (1 initial location + 3 locations to forage an observation from)

d{1} = [1 1 1]';
d{2} = [1 0 0 0]';

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(d);
for f = 1:Nf
    Ns(f) = numel(d{f});
end

No = [2 4]; % two 'what' outcomes (a 'null' outcome and a 'good' outcome), and 4 where outcomes ('where am I foraging from?')
Ng    = numel(No);
for g = 1:Ng
    A{g} = zeros([No(g),Ns]);
end

A{1}(:,:,1) = [1 1 1; 
               0 0 0]; % in the first location, you're not getting anything from any of the hidden states, cuz you're at the start point
A{1}(:,:,2) = [0 0.25 0.75; 
               1 0.75 0.25];
A{1}(:,:,3) = [0.25 0.25 0.25; 
               0.75 0.75 0.75];
A{1}(:,:,4) = [0.9 0.05 0.05;
               0.1 0.95 0.95];
           
A{2}(:,:,1) = [1 1 1;
               0 0 0;
               0 0 0;
               0 0 0;];
A{2}(:,:,2) = [0 0 0;
               1 1 1;
               0 0 0;
               0 0 0;];           
A{2}(:,:,3) = [0 0 0;
               0 0 0;
               1 1 1;
               0 0 0;];           
A{2}(:,:,4) = [0 0 0;
               0 0 0;
               0 0 0;
               1 1 1;]; 
           
           
% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

 
% controllable forage locations: move to the k-th location
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
T         = 6;
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{1}(2,:) =  2;                 % agent wants to encounter observation 2 


% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.T = T;                      % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = d;                      % prior over initial states
mdp.s = [1 1]';                 % initial state indices ([1 1] means state is Goo, and you start in forage location)
mdp.o = [1 1]';                 % initial outcome

mdp.Aname = {'what','where'};
mdp.Bname = {'what','where'};
mdp.temp  = 2;
mdp.alpha = 128;
 
% illustrate a single trial
%==========================================================================
MDP   = spm_MDP_VB_X(mdp);

           


