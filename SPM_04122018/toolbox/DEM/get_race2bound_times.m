function sub_matrix = get_race2bound_times(data_dir,all_files,f_i)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

load(fullfile(data_dir,all_files{f_i}));
num_trials = size(MDPresult,2);


sub_matrix = NaN(num_trials,10);

sub_matrix(:,1) = f_i;


for trial_i = 1:num_trials
    
    sub_matrix(trial_i,2) = trial_i;
    
    current_trial = MDPresult(1,trial_i);
    
    saccades2inspect = find(current_trial.o(1,:) > 1);
    
    time2bound = zeros(1,length(saccades2inspect));
    for s_i = 1:length(saccades2inspect)
        time2bound(s_i) = size(current_trial.mdp(saccades2inspect(s_i)).X{1},2);
    end
    
    sub_matrix(trial_i,saccades2inspect+2) = time2bound;
     
end




end

