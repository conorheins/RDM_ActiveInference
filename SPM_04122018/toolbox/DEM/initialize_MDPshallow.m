function [ MDPshallow ] = initialize_MDPshallow(p)
%INITIALIZE_MDPSHALLOW Initialize lower/shallower level of MDP, with an
%optional manipulation of sensory uncertainty through the p parameter
% INPUTS:  p:           'uncertainty' parameter, used to parameterize the
%                       width  of a Gaussian distribution that introduces 
%                       uncertainty into the mapping between true RDP motion vectors 
%                       and the observed motion (an analogue of RDP motion incoherence)
%                        
% OUTPUTS: MDP_shallow: initialized lower/shallower level of hierarchical MDP

if ~exist('p','var') || isempty(p) || nargin < 1 || p == 0
    cert_vals = [1 0 0 0 ];
else
    cert_vals = normpdf(0:3,0,p);
    cert_vals = [cert_vals(1), cert_vals(2), cert_vals(2), cert_vals(3)];
    cert_vals = cert_vals./sum(cert_vals);
end

D{1} = [ 1 1 1 1 1 ]'; % Direction: {'null','UP','RIGHT','DOWN','LEFT'}

Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end

for f1 = 1:Ns(1)
    if f1 == 1
        A{1}(f1,f1) = 1;
    elseif f1 == 2 % 'UP'
        A{1}(2,f1) = cert_vals(1); % 'UP'
        A{1}(3,f1) = cert_vals(2); % 'RIGHT'
        A{1}(5,f1) = cert_vals(3); % 'LEFT'
        A{1}(4,f1) = cert_vals(4); % 'DOWN'
    elseif f1 == 3 % 'RIGHT'
        A{1}(3,f1) = cert_vals(1); % 'RIGHT'
        A{1}(2,f1) = cert_vals(2); % 'UP'
        A{1}(4,f1) = cert_vals(3); % 'DOWN'
        A{1}(5,f1) = cert_vals(4); % 'LEFT'
    elseif f1 == 4 % 'DOWN'
        A{1}(4,f1) = cert_vals(1); % 'DOWN'
        A{1}(3,f1) = cert_vals(2); % 'RIGHT'
        A{1}(5,f1) = cert_vals(3); % 'LEFT'
        A{1}(2,f1) = cert_vals(4); % 'UP'
    elseif f1 == 5 % 'LEFT'
        A{1}(5,f1) = cert_vals(1); % 'LEFT'
        A{1}(2,f1) = cert_vals(2); % 'UP'
        A{1}(4,f1) = cert_vals(3); % 'DOWN'
        A{1}(3,f1) = cert_vals(4); % 'RIGHT'
    end 
end

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% MDP Structure
%--------------------------------------------------------------------------
MDPshallow.T = 1;                      % number of updates
MDPshallow.A = A;                      % observation model
MDPshallow.B = B;                      % transition probabilities
MDPshallow.D = D;                      % prior over initial states

MDPshallow.Aname = {'Direction'};
MDPshallow.Bname = {'Direction'};
clear A B D
 
MDPshallow = spm_MDP_check(MDPshallow);

end

