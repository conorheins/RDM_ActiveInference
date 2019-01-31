
load RDPsearch_Hierarch15_stats_array.mat

all_sp = [];
all_pp = [];
all_scene_b = {};
all_scene_t = [];

all_sacc_idx = [];

all_EV = [];
all_IV = [];
all_P  = [];

for f_i = 1:length(stats_array)
    
    tic
    temp = stats_array{f_i};
    
    n = size(temp{1},1); % number of trials
            
    for t_i = 1:n
        
        dat = temp{5}{t_i};
        all_EV = [all_EV; ...
            squeeze(dat(:,2:end,1))];
        
        all_IV = [all_IV; ...
            squeeze(dat(:,2:end,2))];
        
        all_P = [all_P; ...
            squeeze(dat(:,2:end,3))];
        
        all_sacc_idx = [all_sacc_idx ; dat(:,1,1)];
        
        all_sp = [all_sp ; ...
            repmat(temp{1}(t_i),size(dat,1),1)];
        
        all_pp = [all_pp ; ...
            repmat(temp{2}(t_i),size(dat,1),1)];
        
        all_scene_b = [all_scene_b ; ...
            repmat(temp{3}(t_i),size(dat,1),1)];
        
        all_scene_t = [all_scene_t ; ...
            repmat(temp{4}(t_i),size(dat,1),1)];
  
    end
    
    fprintf('Time taken to process one file: %.2f seconds\n',toc)
    
end



        
        
        
        
            
    
    
    
    
    
   
        
        
        
        
    
    
    
    
    
    