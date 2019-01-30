function [ MDPshallow ] = initialize_MDPshallow_v6(p,T,policy_depth)
%INITIALIZE_MDPSHALLOW Initialize lower/shallower level of MDP, with an
%optional manipulation of sensory uncertainty through the p parameter
% INPUTS:  p:           'uncertainty' parameter, used to parameterize the
%                       sensory uncertainty into the mapping between true RDP motion vectors 
%                       and the observed motion (an analogue of RDP motion
%                       incoherence). Precisely, p here represents the inverse temperature
%                       parameter of a softmax distribution over perceived
%                       dot directions. The higher p is (the lower the
%                       temperature), the more certain is the mapping
%                       between true dot directions and perceived dot
%                       directions (likelihood matrix approaches identity
%                       matrix with higher p).
%          T:           Number of temporal updates at the level of the MDP, default is 1 
%          policy_depth Temporal horizon of allowable policies (becomes the first
%                       dimension of MDP.V)
%                        
% OUTPUTS: MDP_shallow: initialized lower/shallower level of hierarchical MDP

if ~exist('p','var') || isempty(p) || nargin < 1
    p = 16; % default to high certainty in mapping
end

if ~exist('T','var') || isempty(T) || nargin < 2
    T = 1; % when not provided as argument, T is defaulted to 1 timestep;
end

D{1} = [ 1 1 1 1 1 ]';   % Direction: {'null','UP','RIGHT','DOWN','LEFT'}
D{2} = [ 1 0 0 0 0 0]' ; % Where/choice-state: {'No choice yet','Choose null','Choose UP','Choose RIGHT','Choose DOWN','Choose Left'} 

Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end

% outcome dimensions: 
% 1) visual motion sample: ('null','UP','RIGHT','DOWN','LEFT')
% 2) feedback: {'null','Correct','Incorrect'}
% 3) choice-state: {'No choice yet','Choose null','Choose UP','Choose RIGHT','Choose DOWN','Choose Left'} 

No = [5 3 6];

Ng = numel(No);
for g = 1:Ng
    a{g} = zeros([No(g),Ns]);
end

for f1 = 1:Ns(1)
    
    for f2 = 1:Ns(2)  
        
        if f2 == 1
            
            if f1 == 1
                a{1}(f1,f1,f2) = power(512,10);
            elseif f1 == 2 % 'UP'
                a{1}(2:5,f1,f2) = power(512,10)*spm_softmax([1,0,0,0]'*p);
            elseif f1 == 3 % 'RIGHT'
                a{1}(2:5,f1,f2) = power(512,10)*spm_softmax([0,1,0,0]'*p);
            elseif f1 == 4 % 'DOWN'
                a{1}(2:5,f1,f2) = power(512,10)*spm_softmax([0,0,1,0]'*p);
            elseif f1 == 5 % 'LEFT'
                a{1}(2:5,f1,f2) = power(512,10)*spm_softmax([0,0,0,1]'*p);
            end
            a{2}(1,f1,f2) = power(512,10)*1;
            
        elseif f2 > 1   
            if f1 + 1 == f2
                a{2}(2,f1,f2) = power(512,10)*1; % if your choice (f2) is consistent with the true dot pattern (f1), you're correct
            else
                a{2}(3,f1,f2) = power(512,10)*1; % else, you're wrong
            end
            a{1}(1,f1,f2) = power(512,10)*1; % when you're making your choice, your visual outcome is 'null'
            
            %%%%% DOUBLE CHECK ON THIS -- WILL THIS CAUSE THE AGENT TO
            %%%%% REVISE THEIR BELIEFS & THINK TRUE HIDDEN STATE IS NULL?
            %%%%%% SINCE WHEN THEY'RE CHOOSING THEIR VISUAL OUTCOME IS
            %%%%%% NULL?
        end
        
        % where: choice-state: {'No choice yet','Choose null','Choose UP','Choose RIGHT','Choose DOWN','Choose Left'} 
        %----------------------------------------------------------
        a{3}(f2,f1,f2) = power(512,10)*1;
        
        % this is what it was before
%         a{3}(f2,f1,f2) = 1;
        
    end
end
for g = 1:Ng
    a{g} = double(a{g});
end
A = a; % make the generative process the same as the generative model

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% controllable fixation points: move to the k-th decision
%--------------------------------------------------------------------------
% for k = 1:Ns(2)
%     B{2}(:,:,k) = 0;
%     B{2}(k,:,k) = 1;
% end

% new prior over policies -- decisions ('Choose RIGHT',etc.) are 'sinks' in
% the transition matrix. I.e. once you decide, there's no going back
% (21.01.2019)
%--------------------------------------------------------------------------
for k = 2:Ns(2)
    B{2}(:,2:end,k) = B{2}(:,2:end,1);
    B{2}(k,1,k) = 1;
end

% allowable policies (here, specified as the next action) U
%--------------------------------------------------------------------------

if ~exist('policy_depth','var') || isempty(p) || nargin < 3
    policy_depth = 1; % default to one-step ahead policies
end

% allowable policies (here, specified as the next action) U
%--------------------------------------------------------------------------
Np        = size(B{2},3);
U         = ones(policy_depth,Np,Nf);
U(:,:,2)  = 1:Np;

% priors: (utility) C
%--------------------------------------------------------------------------
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{3}      = zeros(No(3),T);
% C{2}(1,:) =  2;               % option to add a cost to being undecided
C{2}(2,:) =  2;                 % the agent expects to be right
C{2}(3,:) = -4;                 % and not wrong

% MDP Structure
%--------------------------------------------------------------------------
MDPshallow.T = T;                      % number of updates
MDPshallow.U = U;                      % allowable policies
MDPshallow.a = a;                      % concentration parameters for observation model
MDPshallow.A = A;                      % observation model
MDPshallow.B = B;                      % transition probabilities
MDPshallow.C = C;                      % preferred outcomes
MDPshallow.D = D;                      % prior over initial states

MDPshallow.Aname = {'Direction','Feedback','Choice-State'};
MDPshallow.Bname = {'Direction','Choice-State'};
 
MDPshallow = spm_MDP_check(MDPshallow);

end

