function [ MDPshallow ] = initialize_trustMDP_shallow(k,U)
%INITIALIZE_TRUSTMDP_SHALLOW Initialize lower/shallower level of MDP
% INPUTS:  k:           edge parameter that determines number of players
%                       that you will play with
%          U:           cell array of payoff matrices 
%                        
% OUTPUTS: MDP_shallow: initialized lower/shallower level of hierarchical trust MDP

% initial state - encoding the start position and the two possible actions
% of your partner
%--------------------------------------------------------------------------
% S    = [1 0 0 0 0]';        % true state of the world (to be used in the generative process)

% prior beliefs about initial state
%--------------------------------------------------------------------------
k    = [1 1]; % whether the game partner is social/anti-social
p    = spm_softmax(k(:));
k    = log(p);
D    = kron(p,[1 0 0 0 0]');

% investor's payoffs or prior beliefs (softmax(utility))
%--------------------------------------------------------------------------
a    = 1/2; % kinda encodes 1) your altruism and 2) your beliefs about {altruism of your opponent |  they're pro-social}
pp   = [0; spm_softmax(spm_vec((1 - a)*U{1} + a*U{2}))]; % encodes preferences, since higher values of U{1} -- YOUR payoff matrix --
% get higher probability. High payoffs of your opponent/trustee also get
% high probability if you're altruism (a) is high/close to 1.

pn   = [0; spm_softmax(spm_vec(        U{1}         ))];

C    = [pp*p(1); pn*p(2)];

% investor's belief (based on a prosocial and nonsocial trustee)
%--------------------------------------------------------------------------
cp   = spm_softmax(((1 - a)*U{2}(1,:) + a*U{1}(1,:))'); % your beliefs about {possibilities of CC/CD | they're pro-social}
dp   = spm_softmax(((1 - a)*U{2}(2,:) + a*U{1}(2,:))'); % your beliefs about {possibilities of DC/DD | they're pro-social}
cn   = spm_softmax((        U{2}(1,:)              )'); % your beliefs about {possibilities of CC/CD | they're anti-social}
dn   = spm_softmax((        U{2}(2,:)              )'); % your beliefs about {possibilities of DC/DD | they're anti-social}



MDPshallow = spm_MDP_check(MDPshallow);

end