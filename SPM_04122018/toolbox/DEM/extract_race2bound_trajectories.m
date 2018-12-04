%% extract_race2bound trajectories (cumulative free energy)

% choose data with MDP structure to analyze

data_dir = uigetdir();

all_files = dir(fullfile(data_dir,'*.mat'));
all_files = {all_files(:).name};

num_files = length(all_files);

master_array = {};
master_metadata = [];

tic
for f_i = 1:num_files
    
    [sub_array,sub_metadata] = get_race2bound_curves(data_dir,all_files,f_i);
       
    master_array = [master_array;sub_array];
    
    master_metadata = [master_metadata;sub_metadata];
    
    clear sub_array sub_metadata;
    
end
fprintf('Time taken to load and extract race2bound curves from all %d files: %.2f minutes\n',num_files,toc/60)

%% gather evidence-accumulation trajectories (i.e. time-dependent accumulation of negative free energy), for different sensory precisions/prior precisions

coherences = master_metadata(:,3);
prior_precisions = master_metadata(:,4);

% coherence_bins = exp(linspace(log(min(coherences)),log(max(coherences)),11));
% coherence_bins = linspace(min(coherences),max(coherences),20);
coherence_bins = linspace(min(coherences),4,11);

% coherence_bins = exp(linspace(log(min(coherences)),log(5),9));
% p_bins = exp(linspace(log(min(prior_precisions)),log(max(prior_precisions)),2));
p_bins = [0.25 0.55 0.8];
% p_bins = linspace(min(prior_precisions),max(prior_precisions),7);

all_trajectories = zeros(length(coherence_bins)-1,length(p_bins)-1,32,2);

for coh_i = 1:length(coherence_bins)-1
    
    coh_bin_idx = and(master_metadata(:,3) >= coherence_bins(coh_i),master_metadata(:,3) < coherence_bins(coh_i+1));
    
    for p_i = 1:length(p_bins)-1
        
        p_bin_idx = and(master_metadata(:,4) >= p_bins(p_i),master_metadata(:,4) < p_bins(p_i+1));
        
        intersection = coh_bin_idx&p_bin_idx;
        
        trials2inspect = master_array(intersection);
        
        trajs = [];
        
        for traj_i = 1:length(trials2inspect)
            if size(trials2inspect{traj_i}{1,1},2) == 32
%                 trajs = [trajs;trials2inspect{traj_i}{1,1}];
                trajs = [trajs;-trials2inspect{traj_i}{2,1}];
            end
        end
        
        if ~isempty(trajs)
            all_trajectories(coh_i,p_i,:,1) = mean(trajs,1);
            all_trajectories(coh_i,p_i,:,2) = std(trajs,0,1)./size(trajs,1);
%             all_trajectories(coh_i,p_i,:,2) = std(trajs,0,1);
        end
        
    end
end


%% plotting

vectors = spm_softmax([coherence_bins;zeros(3,length(coherence_bins))]);
coherence_axis = vectors(1,:);          

coherence_labels = cell(1,length(coherence_axis)-1);
for i = 1:length(coherence_axis)-1
    coherence_labels{1,i} = sprintf('Coherence: %.1f to %.1f %%',coherence_axis(i)*100,coherence_axis(i+1)*100);
end

coherence_colors = cool(length(coherence_labels));
coherence_colors = flipud(coherence_colors(1:length(coherence_labels),:));
for i = 1:size(coherence_colors,1)
    lineProps.col{i} = coherence_colors(i,:);
end


figure;
for p_i = 1:length(p_bins)-1
    
    mseb(1:32,squeeze(all_trajectories(:,p_i,:,1)),squeeze(all_trajectories(:,p_i,:,2)),lineProps);
%     mseb(3:10,squeeze(all_trajectories(:,p_i,3:10,1)),squeeze(all_trajectories(:,p_i,3:10,2)),lineProps);
%     mseb(1:6,squeeze(all_trajectories(:,p_i,1:6,1)),squeeze(all_trajectories(:,p_i,1:6,2)),lineProps);


    legend(coherence_labels);
    
    title(sprintf('Evidence accumulation curves for different coherences, prior precisions between %.2f and %.2f %%',p_bins(p_i)*100,p_bins(p_i+1)*100),...
        'FontSize',24);
    
    xlim([1 32])
    xlabel('Temporal sample index','FontSize',20)
    ylabel('Negative free energy','FontSize',20)
    xt = get(gca,'XTickLabel');
    yt = get(gca,'YTickLabel');
    set(gca,'XTickLabel',xt,'fontsize',16);
    set(gca,'YTickLabel',yt,'fontsize',16);
    
    pause;
    

    close gcf;
    
end

%% look at free energy curves as a function of time

saccades2inspect = 2:8;

columns = 5 + saccades2inspect; % 5 columns in master_metadata before saccade columns

coherences = master_metadata(:,3);
coh = master_metadata(:,3); % get sensory precisions (to be converted to coherences)
coh_bins = linspace(min(coh),3,5);

vectors = spm_softmax([coh_bins;zeros(3,length(coh_bins))]);
coh_axis = vectors(1,:);          

coh_labels = cell(1,length(coh_axis)-1);
for i = 1:length(coh_axis)-1
    coh_labels{1,i} = sprintf('Coherence: %.1f to %.1f %%',coh_axis(i)*100,coh_axis(i+1)*100);
end

coh_colors = cool(length(coh_labels));
coh_colors = flipud(coh_colors(1:length(coh_labels),:));
for i = 1:size(coh_colors,1)
    lineProps.col{i} = coh_colors(i,:);
end

for sacc_i = 1:length(saccades2inspect)
    
    CoI = columns(sacc_i);
    
    all_trajs = zeros(length(coh_bins)-1,32,2);
    for coh_i = 1:length(coh_bins)-1
        
        intersection = ~isnan(master_metadata(:,CoI))&and(coherences >= coh_bins(coh_i),coherences<coh_bins(coh_i+1))...
            &master_metadata(:,CoI)==32;
        
        if any(istrue(intersection))
            trials2inspect = master_array(intersection);
            
            trajs = [];
            
            for traj_i = 1:length(trials2inspect)
                
                sacc_idx = cell2mat(trials2inspect{traj_i}(3,:));
                trajs = [trajs;trials2inspect{traj_i}{1,find(sacc_idx == saccades2inspect(sacc_i))}];
                %             trajs = [trajs;-trials2inspect{traj_i}{2,1}];
            end
            
            all_trajs(coh_i,:,1) = mean(trajs,1);
            all_trajs(coh_i,:,2) = std(trajs,0,1)./sqrt(size(trajs,1));
        end
        
    end
    
    figure(sacc_i);
    mseb(1:32,squeeze(all_trajs(:,:,1)),squeeze(all_trajs(:,:,2)),lineProps);
    title(sprintf('Evidence accumulation for difference coherences, saccade %d',saccades2inspect(sacc_i)))
    
    pause; close gcf;
    
end

    

   
    
    
    
    
    
    
    
    
    



