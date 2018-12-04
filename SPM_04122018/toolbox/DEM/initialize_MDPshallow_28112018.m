function [ MDPshallow ] = initialize_MDPshallow_28112018(p,T,policy_depth)
%INITIALIZE_MDPSHALLOW Initialize lower/shallower level of MDP, with an
%optional manipulation of sensory uncertainty through the p parameter
% INPUTS:  p:           'uncertainty' parameter, used to parameterize the
%                       width  of a Gaussian distribution that introduces 
%                       uncertainty into the mapping between true RDP motion vectors 
%                       and the observed motion (an analogue of RDP motion incoherence)
%          T:           Number of temporal updates at the level of the MDP, default is 1 
%          policy_depth Temporal horizon of allowable policies (becomes the first
%                       dimension of MDP.V)
%                        
% OUTPUTS: MDP_shallow: initialized lower/shallower level of hierarchical MDP

if ~exist('p','var') || isempty(p) || nargin < 1
    p = 16; % default to total certainty in mapping
end

if ~exist('T','var') || isempty(T) || nargin < 2
    T = 1; % when not provided as argument, T is defaulted to 1 timestep;
end

D{1} = [ 1 1 1 1 1 ]'; % Direction: {'null','UP','RIGHT','DOWN','LEFT'}

Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end

for f1 = 1:Ns(1)
    if f1 == 1
        a{1}(f1,f1) = power(512,10);       
    elseif f1 == 2 % 'UP'
        a{1}(2:5,f1) = power(512,10)*spm_softmax([1,0,0,0]'*p);       
    elseif f1 == 3 % 'RIGHT'
        a{1}(2:5,f1) = power(512,10)*spm_softmax([0,1,0,0]'*p);
    elseif f1 == 4 % 'DOWN'
        a{1}(2:5,f1) = power(512,10)*spm_softmax([0,0,1,0]'*p);
    elseif f1 == 5 % 'LEFT'
        a{1}(2:5,f1) = power(512,10)*spm_softmax([0,0,0,1]'*p);
    end 
end

A = a; % make the generative process the same as the generative model

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

 
% allowable policies (here, specified as the next action) U
%--------------------------------------------------------------------------

if ~exist('policy_depth','var') || isempty(p) || nargin < 3
    policy_depth = 1; % default to total certainty in mapping
end

Np        = size(B{1},3);
U         = ones(policy_depth,Np,Nf);
U(:,:,1)  = 1:Np;

% MDP Structure
%--------------------------------------------------------------------------
MDPshallow.T = T;                      % number of updates
MDPshallow.U = U;                      % allowable policies
MDPshallow.a = a;                      % concentration parameters for observation model
MDPshallow.A = A;                      % observation model
MDPshallow.B = B;                      % transition probabilities
MDPshallow.D = D;                      % prior over initial states

MDPshallow.Aname = {'Direction'};
MDPshallow.Bname = {'Direction'};
 
MDPshallow = spm_MDP_check(MDPshallow);

end

