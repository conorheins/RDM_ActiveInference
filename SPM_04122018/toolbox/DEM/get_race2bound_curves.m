function [sub_array, sub_metadata] = get_race2bound_curves(data_dir,all_files,f_i)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

load(fullfile(data_dir,all_files{f_i}));
num_trials = size(MDPresult,2);

sub_array = cell(num_trials,1);

sub_metadata = NaN(num_trials,13);
sub_metadata(:,1) = f_i;
sub_metadata(:,2) = 1:num_trials;
sub_metadata(:,3:4) = [noise_vals,prior_scene_probs];

for trial_i = 1:num_trials
   
    current_trial = MDPresult(1,trial_i);
    
    sub_metadata(trial_i,5) = true_scenes(trial_i) == current_trial.s(1,1);

    saccades2inspect = find(current_trial.o(1,:) > 1);
    
    time2bound = zeros(1,length(saccades2inspect));
    
    sub_array{trial_i,1} = cell(2,length(saccades2inspect));
    for s_i = 1:length(saccades2inspect)
        time2bound(s_i) = size(current_trial.mdp(saccades2inspect(s_i)).X{1},2);
        sub_array{trial_i,1}{1,s_i} = current_trial.mdp(saccades2inspect(s_i)).H_int; % get the cumulative free energy timecourse
        sub_array{trial_i,1}{2,s_i} = current_trial.mdp(saccades2inspect(s_i)).H; % get the free energy timecourse
        sub_array{trial_i,1}{3,s_i} = saccades2inspect(s_i); % record the time-index of the saccade
    end
    
    sub_metadata(trial_i,saccades2inspect+5) = time2bound;
     
end


end

