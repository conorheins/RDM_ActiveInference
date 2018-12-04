%% collect_stats_across_trials function that is called within master loop
function sub_matrix = collect_stats_across_trials(data_dir,all_files,file_idx)

load(fullfile(data_dir,all_files{file_idx}));
num_trials = size(MDPresult,2);
sub_matrix = zeros(num_trials,6);

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
        sub_matrix(trial_i,5) = length(most_believed);
    else
        sub_matrix(trial_i,5) = find(most_believed == true_state_idx,1);
    end
    
    % calculate repeated visits to quadrants
    o = current_trial.o(3,:);
    sub_matrix(trial_i,6) = length(o(o<6)) - length(unique(o(o<6)));
end



end

