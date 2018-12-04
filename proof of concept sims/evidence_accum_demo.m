function [means,sems] = evidence_accum_demo(num_outcomes,num_samples,num_coherences,num_iter)
% EVIDENCE_ACCUM_DEMO:
% Inputs: num_outcomes (e.g. 4 possible sensory outcomes); num_samples
% (e.g. 32 independent samples from a probability distribution); a
% number of coherences/sensory precisions to iterate over (num_coherences);
% and number of iterations to get statistics over the belief-updating for each coherence level (num_iter)

coherence_values = linspace(1/num_outcomes,1,num_coherences); % set of coherence values with which to parameterize sensory uncertainty 

% plotting stuff
coherence_colors = cool(ceil(1.5*length(coherence_values)));
coherence_colors = flipud(coherence_colors(1:length(coherence_values),:));

%initialize likelihood matrix
init_likelihood = eye(num_outcomes); % initial mapping from hidden states (true RDM motion vectors) to outcomes (perceived/sampled RDM motion vectors)


means = zeros(length(coherence_values),num_samples);
sems = zeros(length(coherence_values),num_samples);

figure;
for i = 1:length(coherence_values)
    
    coherence = coherence_values(i);
    
    likelihood = init_likelihood;
    
    % fill out likelihood matrix according to coherence
    for j = 1:num_outcomes
        current_column = logical(likelihood(:,j));
        likelihood(current_column,j) = coherence;
        likelihood(~current_column,j) = (1 - coherence)/(num_outcomes-1);    
    end
    
    % randomly sample a true state
    [~,true_state] = max(rand(num_outcomes,1)); % index of the true state (1 == 'UP', 2 == 'RIGHT', 3 == 'DOWN', 4 == 'LEFT', etc.)
   
    log_likelihood = log(likelihood + exp(-16)); % add small values to avoid Inf in the log
   
    all_trajectories = zeros(num_iter,num_samples);
    
    for iter = 1:num_iter
        % generate sensory data by sampling from the appropriate column of the
        % likelihood mapping
        obs = zeros(1,num_samples);
        for j = 1:length(obs)
            obs(j) = find(rand < cumsum(likelihood(:,true_state)),1);
        end
        
        % start with flat prior beliefs over dot directions
        prior_beliefs = softmax(ones(num_outcomes,1));
        
        % initialize posterior beliefs to prior beliefs
        posterior_beliefs = prior_beliefs;
        
        belief_history = zeros(num_outcomes,num_samples);
        
        for sample = 1:num_samples
            posterior_beliefs = log_likelihood(obs(sample),:)' + log(posterior_beliefs); % Bayesian belief update
            posterior_beliefs = softmax(posterior_beliefs); % re-normalize to make a probability distribution
            belief_history(:,sample) = posterior_beliefs; % store history of updates
        end
        
        all_trajectories(iter,:) = belief_history(true_state,:);
        
    end
    
    means(i,:) = mean(all_trajectories,1);
    sems(i,:) = std(all_trajectories,0,1)./sqrt(num_iter);
    
end

coh_names = cell(1,length(coherence_values));
for coh_i = 1:length(coherence_values)
    lineProps.col{coh_i} = coherence_colors(coh_i,:);
    coh_names{coh_i} = sprintf('Coherence value: %.1f %%',coherence_values(coh_i)*100);
end

mseb(1:num_samples,means,sems,lineProps);


xlim([1 num_samples]);
xlabel(sprintf('Sequential sample index (%d total)',num_samples))
ylabel('Posterior belief in true state')
legend(coh_names)
legend('show')

function [y] = softmax(x,k)
% softmaxfunction of COLUMN vectors
% FORMAT [y] = spm_softmax(x,k)
%
% x - vector of activity
% k - temperature or inverse sensitivity parameter (default k = 1)
%
% y   = exp(k*x)/sum(exp(k*x))
%
% NB: If supplied with a matrix this routine will return the softmax
% function over columns - so that spm_softmax([x1,x2,..]) = [1,1,...]
 
% apply
%--------------------------------------------------------------------------
if nargin > 1, x = k*x; end

n    = size(x);
if n(end) > 1
    ind    = cell(size(n));
    ind(:) = {':'};
    for num_i  = 1:n(end)
        sub       = ind;
        sub(end)  = {num_i};
        y(sub{:}) = softmax(x(sub{:}));
    end
else
    x   = x - max(x);
    ex  = exp(x);
    y   = ex/sum(ex);
end

end

end

