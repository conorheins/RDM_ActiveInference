%% solve and evaluate performance for sequence of trials

N = 100; % number of trials
s(1,:) = ceil(rand(1,N)*48);            % randomly initialize states
noise_val(1,:) = 4.*rand(1,N);          % randomly initiate coherence parameter
scene_idx = ceil(rand(1,N) * 4); 
prior_scenes = cell(1,N);               % which scene category does the agent believe occupies most trials, a priori
prior_scenes(scene_idx == 1) = {'UP_RIGHT'};
prior_scenes(scene_idx == 2) = {'RIGHT_DOWN'};
prior_scenes(scene_idx == 3) = {'DOWN_LEFT'};
prior_scenes(scene_idx == 4) = {'LEFT_UP'};
prior_scene_prob = 0.25 + (1-0.25).*rand(1,N);                  % agent's prior belief in proportion of trials that have scene identity 'prior_scene'

for i = 1:N
    [MDP(i),Scenes]   = initialize_searchMDP_4(noise_val(i),s(i),prior_scenes{i},prior_scene_prob(i));      % create structure array
end

tic
MDP = spm_MDP_VB_X_distribO(MDP);
fprintf('Time taken to run a %d trials: %.2f minutes\n',N,toc/60);

for i = 1:N
    view_trial_trajectory(MDP,i,Scenes);
    fprintf('Trial %d, uncertainty value: %.2f, prior certainty in scene %s: %.2f\n',i,noise_val(i),prior_scenes{i},prior_scene_prob(i));
    pause; 
    close;
end

for i = 1:N
    o      = MDP(i).o(1,:);                       % outcomes
    stats_matrix(1,i) = double(any(o == 6) & ~any(o == 7));  % accuracy
    stats_matrix(2,i) = find([(o > 5), 1],1) - 1;            % number of saccades
    stats_matrix(3,i) = mean(MDP(i).rt);                     % reaction time
end

[counts_incor,bins_incor] = histcounts(noise_val(stats_matrix(1,:)==0),10,'Normalization','probability');
[counts_cor,bins_cor] = histcounts(noise_val(stats_matrix(1,:)==1),10,'Normalization','probability');


save(sprintf('MDPresults_%s.mat',datestr(now,'mmddyy')),'MDP','-v7.3')