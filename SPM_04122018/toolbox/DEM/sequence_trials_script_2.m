%% solve and evaluate performance for sequence of trials

N = 10; % number of trials
s(1,:) = ceil(rand(1,N)*48); % randomly initialize states
p(1,:) = rand(1,N);          % randomly initiate coherence parameter

for i = 1:N
    [MDP(i),Scenes]   = initialize_searchMDP_2(p(i),s(i));      % create structure array
end

tic
MDP = spm_MDP_VB_X(MDP);
fprintf('Time taken to run a %d trials: %.2f minutes\n',N,toc/60);

for i = 1:N
    view_trial_trajectory(MDP,i,Scenes);
    fprintf('Trial %d, uncertainty value: %.2f\n',i,p(i));
    pause; 
    close;
end

for i = 1:N
    o      = MDP(i).o(1,:);                       % outcomes
    stats_matrix(1,i) = double(any(o == 6) & ~any(o == 7));  % accuracy
    stats_matrix(2,i) = find([(o > 5), 1],1) - 1;            % number of saccades
    stats_matrix(3,i) = mean(MDP(i).rt);                     % reaction time
end

[counts_incor,bins_incor] = histcounts(p(stats_matrix(1,:)==0),10,'Normalization','probability');
[counts_cor,bins_cor] = histcounts(p(stats_matrix(1,:)==1),10,'Normalization','probability');


save(sprintf('MDPresults_%s.mat',datestr(now,'mmddyy')),'MDP','-v7.3')