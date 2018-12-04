%% collect_stats_across_trials function that is called within master loop
function sub_matrix = collect_stats_across_trialsv2(data_dir,all_files,file_idx)

load(fullfile(data_dir,all_files{file_idx}));
num_trials = size(MDPresult,2);

% statistics to collect:
% 1. noise value (width of the gaussian that is filling out the columns of
%    the likelihood matrix mapping from dot patterns to incoherent realizations thereof -- proxy for sensory uncertainty)
% 2. prior confidence of beliefs in 1 of 4 scenes, believed by the agent
%    at the outset of the trial
% 3. true state (index from 1 to 4)
% 4. initial believed state (index from 1 to 4)
% 5. time until 'realization' of true state
% 6. number of repeated visits to quadrants (not including revisits that
%    happen after categorization
% 7. Accuracy ( 1 = Correct, 2 = Incorrect, 3 = Undecided)
% 8. Biggest prediction error value (difference in free energies between
% 	 previous and current time, biggest moment at which that happened)
% 9. Time of biggest prediction error value
% 10. Average entropy of likelihood over outcomes (before categorization
%    decision, and only those quadrants with a dot pattern in them)
% 11. When did was the highest entropy of likelihood over outcomes
%   encountered
% 12. When did was the lowest entropy of likelihood over outcomes
%   encountered


sub_matrix = zeros(num_trials,12);

sub_matrix(:,1) = noise_vals;
sub_matrix(:,2) = prior_scene_probs;

for trial_i = 1:num_trials
    
    current_trial = MDPresult(1,trial_i);
    
    posteriors_about_states = zeros(size(current_trial.xn{1},2),size(current_trial.xn{1},3));
    for timestep = 1:size(posteriors_about_states,2)
        posteriors_about_states(:,timestep) = squeeze(current_trial.xn{1}(end,:,timestep,timestep));
    end
    
    prior_strengths = zeros(4,1);
    for scene_i = 1:4
        if ismember(current_trial.s(1,1),scene_idx{scene_i})
            true_state_idx = scene_i;
        end
        prior_strengths(scene_i) = sum(current_trial.D{1}(scene_idx{scene_i}));
    end
    
    sub_matrix(trial_i,3) = true_state_idx;
    [~,sub_matrix(trial_i,4)] = max(prior_strengths);
    
    most_believed = zeros(1,size(posteriors_about_states,2));
    for timestep = 1:size(posteriors_about_states,2)
        
        each_scene_beliefs = zeros(4,1);
        beliefs_at_t = posteriors_about_states(:,timestep);
        for scene_i = 1:4
            each_scene_beliefs(scene_i,1) = sum(beliefs_at_t(scene_idx{scene_i}));
        end
        [~,most_believed(timestep)] = max(each_scene_beliefs);
        
    end
    
    if isempty(find(most_believed == true_state_idx,1))
        sub_matrix(trial_i,5) = NaN; % they never realized it
    else
        sub_matrix(trial_i,5) = find(most_believed == true_state_idx,1);
    end        
    
    % calculate repeated visits to quadrants, only including part of trial
    % before categorization
    o = current_trial.o(3,:);
    pre_cat_time = find(o>=6,1); % timestep of categorization
    if isempty(pre_cat_time) % if they never decided/categorized
        sub_matrix(trial_i,6) = length(o(o<6)) - length(unique(o(o<6)));
    elseif pre_cat_time > 1 % if they took more than one saccade to decide
        o = o(1:pre_cat_time-1);
        sub_matrix(trial_i,6) = length(o) - length(unique(o));
    else % if they decided at the first timestep
        sub_matrix(trial_i,6) = 0;
    end
    
    % second row of o has the accuracy read-out: 1 if correct, 2 if
    % incorrect, 3 if undecided
    o2 = current_trial.o(2,:);
    if any(o2 == 2) && ~any(o2 == 3)
        sub_matrix(trial_i,7) = 1;
    elseif ~any(o2 == 2 | o2 == 3)
        sub_matrix(trial_i,7) = 3;
    else
        sub_matrix(trial_i,7) = 2;
    end
    
 
    o = current_trial.o(3,:);
    pre_cat_time = find(o>=6,1); % timestep of categorization
    if isempty(pre_cat_time) || pre_cat_time == 1 % if they never decided/categorized
        free_energies = current_trial.G(6:end,1:end-1);
    elseif pre_cat_time > 1 % if they took more than one saccade to decide
        free_energies = current_trial.G(6:end,1:pre_cat_time-1);
    end
    
    % track expected free energies of each categorization, calculate
    % the change in free energy 
    best_policy = zeros(1,size(free_energies,2));
    
    for t = 1:size(free_energies,2)
        [~,best_policy(t)] = max(free_energies(:,t));
    end
    
    belief_changes = find(diff(best_policy) ~= 0);
    if isempty(belief_changes)
        if isempty(pre_cat_time)
            free_energy_diffs = abs(diff(current_trial.H));
            [sub_matrix(trial_i,8),best_id] = max(free_energy_diffs);
            sub_matrix(trial_i,9) = best_id + 1;
        elseif pre_cat_time == 2
            sub_matrix(trial_i,8) = abs(diff(current_trial.H(1:2)));
            sub_matrix(trial_i,9) = 2;
        else
            free_energy_diffs = abs(diff(current_trial.H(1:(pre_cat_time-1))));
            [sub_matrix(trial_i,8),best_id] = max(free_energy_diffs);
            sub_matrix(trial_i,9) = best_id + 1;
        end
    else
        belief_changes = belief_changes + 1;
        free_energy_diffs = [];
        for i = 1:length(belief_changes)
            grad = abs(current_trial.H(belief_changes(i)) - current_trial.H(belief_changes(i)-1));
            free_energy_diffs = [free_energy_diffs,grad];
        end
        [sub_matrix(trial_i,8),best_id] = max(free_energy_diffs);
        sub_matrix(trial_i,9) = belief_changes(best_id);
    end
    
    if isempty(pre_cat_time)
        T_indices = find(and(current_trial.o(3,:) > 1, current_trial.o(3,:) < 6));
    else
        T_indices = find(and(current_trial.o(3,1:pre_cat_time-1) > 1, current_trial.o(3,1:pre_cat_time-1) < 6));
    end
    
    if isempty(T_indices)
        sub_matrix(trial_i,10:12) = NaN;
    else 
        all_Ox = zeros(1,length(T_indices));
        for t = 1:length(T_indices)
            LO = current_trial.O{1,T_indices(t)}(2:5);
            all_Ox(t) = -LO'*log(LO + 1e-16);
        end
        
        sub_matrix(trial_i,10) = mean(all_Ox);
        [~,sub_matrix(trial_i,11)] = max(all_Ox);
        [~,sub_matrix(trial_i,12)] = min(all_Ox);
    end

    
    
end
    
    
end    

    
      