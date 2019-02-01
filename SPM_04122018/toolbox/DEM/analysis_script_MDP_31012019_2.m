
load RDPsearch_Hierarch15_stats_array.mat


%% accumulate data from stats array into big matrices (makes easier to index by variable/condition)
all_sp = [];
all_pp = [];
all_scene_b = {};
all_scene_t = [];

all_sacc_idx = [];

all_EV = [];
all_IV = [];
all_P  = [];

time_to_analyze = zeros(1,length(stats_array));

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
    
    time_to_analyze(f_i) = toc;
    
end


%% plot trajectories of characteristic quantities (differential epistemic value over sampling time, instrumental value over sampling time, probability of deciding over sampling time)

sp_bins = linspace(min(all_sp),max(all_sp),6);
pp_bins = linspace(min(all_pp),max(all_pp),6);

for sp_i = 1:length(sp_bins)-1
    
    sp_bin_idx = all_sp >= sp_bins(sp_i) & all_sp < sp_bins(sp_i+1);
    
    time2break_bins = [2 6 12 20];
    
    times2break = zeros(size(all_EV,1),1);
    
    for ii = 1:size(all_EV,1)
        temp = find(isnan(all_EV(ii,:)),1);
        if isempty(temp)
            times2break(ii) = T;
        else
            times2break(ii) = temp;
        end
    end
    
     for tb_i = 1:length(time2break_bins)-1
         
          tb_bin_idx = times2break >= time2break_bins(tb_i) & times2break < time2break_bins(tb_i+1);
               
          for pp_i = 1:length(pp_bins)-1
              
              pp_bin_idx = all_pp >= pp_bins(pp_i) & all_pp < pp_bins(pp_i+1);
              
              EV_matches = all_EV(sp_bin_idx & tb_bin_idx & pp_bin_idx,:);
              
              IV_matches = all_IV(sp_bin_idx & tb_bin_idx & pp_bin_idx,:);
              
              P_matches = all_P(sp_bin_idx & tb_bin_idx & pp_bin_idx,:);
              
              plot(mean(EV_matches,1,'omitnan'),'DisplayName',sprintf('Prior beliefs between %.0f and %.0f %%',...
                  pp_bins(pp_i)*100,pp_bins(pp_i+1)*100));
              
%               plot(mean(IV_matches,1,'omitnan'),'DisplayName',sprintf('Prior beliefs between %.0f and %.0f %%',...
%                   pp_bins(pp_i)*100,pp_bins(pp_i+1)*100));
%               plot(mean(P_matches,1,'omitnan'),'DisplayName',sprintf('Prior beliefs between %.0f and %.0f %%',...
%                   pp_bins(pp_i)*100,pp_bins(pp_i+1)*100));
                         
              
              hold on;
              
          end
          
          legend('show');
          
          pause; close gcf;
        
    end
    
end
        
        
        
        
        
        
    





        
        
        
        
            
    
    
    
    
    
   
        
        
        
        
    
    
    
    
    
    