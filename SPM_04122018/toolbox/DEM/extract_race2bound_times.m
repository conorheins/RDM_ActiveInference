%% extract_race2bound_times

% choose data with MDP structure to analyze

data_dir = uigetdir();

all_files = dir(fullfile(data_dir,'*.mat'));
all_files = {all_files(:).name};

num_files = length(all_files);

indicator_matrix = [];

tic
for f_i = 1:num_files
    
    sub_matrix = get_race2bound_times(data_dir,all_files,f_i);
       
    indicator_matrix = [indicator_matrix;sub_matrix];
    
    clear sub_matrix;
    
end
fprintf('Time taken to load and extract race2bound_times from all %d files: %.2f minutes\n',num_files,toc/60)

%%

list_of_indices = [];

for i = 1:size(indicator_matrix,1)
    if any(and(indicator_matrix(i,3:end) > 4,indicator_matrix(i,3:end) < 32))
        list_of_indices = [list_of_indices;indicator_matrix(i,:)];
    end
end

f_ids = unique(list_of_indices(:,1));

chi = 1/10;

for f_i = 1:length(f_ids)
    
    load(fullfile(data_dir,all_files{f_ids(f_i)}));
    
    trials_to_inspect = list_of_indices(list_of_indices(:,1) == f_ids(f_i),:);
    
    for trial_i = 1:size(trials_to_inspect,1)
        
        saccades2inspect = find(and(trials_to_inspect(trial_i,3:end) > 4, trials_to_inspect(trial_i,3:end)<32));
        
        for s_i = 1:length(saccades2inspect)
            
            H_int = MDPresult(trials_to_inspect(trial_i,2)).mdp(saccades2inspect(s_i)).H_int;
            plot(H_int)
            hold on; plot(1:length(H_int),-chi*ones(1,length(H_int)))
            
            title(sprintf('FileID %d, Trial %d, Saccade %d',f_ids(f_i),trials_to_inspect(trial_i,2),saccades2inspect(s_i)))
            
            pause;
            close gcf;
            
        end
        
    end
    
end
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



