% analysis script 31012019


%% set paths/directories

results_dir = uigetdir();

results_fnam = fullfile(results_dir,uigetfile('*.mat'));


%% load data 

load(results_fnam)

%% analyze

var_names = {'Sensory Precision', 'Prior Belief Precision', 'True Scene','Scene Belief','SaccadeNum','Time to Choose'}; 

sPrecisions = results_array(:,strcmp(var_names,'Sensory Precision'));
pPrecisions = results_array(:,strcmp(var_names,'Prior Belief Precision'));


%% plot latency to choose as a function of sensory and prior precision

% sp_bins = linspace(min(sPrecisions),max(sPrecisions),11);

% sp_bins = linspace(0.5,1.5,6);

% sp_bins = exp(linspace(log(1),log(max(sPrecisions)),11));
% sp_bins(2:4) = [];

sp_bins = linspace(1.7,max(sPrecisions),11);
sp_bins = exp(linspace(log(1),log(6),8));

pp_bins = linspace(min(pPrecisions),max(pPrecisions),7);
% pp_bins = exp(linspace(log(min(pPrecisions)),log(max(pPrecisions)),5));

latency_matrix = zeros(length(sp_bins)-1,length(pp_bins)-1,2);

lat_col_idx = strcmp(var_names,'Time to Choose');

for sp_i = 1:length(sp_bins)-1
    
    sp_bin_idx = and(results_array(:,strcmp(var_names,'Sensory Precision')) >= sp_bins(sp_i),results_array(:,strcmp(var_names,'Sensory Precision')) < sp_bins(sp_i+1)); 
    
    for pp_i = 1:length(pp_bins)-1
        p_bin_idx = and(results_array(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(pp_i),results_array(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(pp_i+1));
        latencies = results_array(sp_bin_idx&p_bin_idx,lat_col_idx); % latencies
        latency_matrix(sp_i,pp_i,1) = mean(latencies);
        latency_matrix(sp_i,pp_i,2) = std(latencies)./sqrt(length(latencies));
%         latency_matrix(sp_i,pp_i,1) = median(latencies);
%         latency_matrix(sp_i,pp_i,2) = mad(latencies);
    end
       
end

latency_means = squeeze(latency_matrix(:,:,1));
latency_errors = squeeze(latency_matrix(:,:,2));
sp_axis = sp_bins(1:end-1) + diff(sp_bins)/2;
prior_precis = 100*(pp_bins(1:end-1) + diff(pp_bins)/2);
prior_precision_labels = cell(1,length(pp_bins)-1);
for i = 1:length(prior_precision_labels)
    prior_precision_labels{1,i} = sprintf('%.0f - %0.f %% prior belief',pp_bins(i)*100,pp_bins(i+1)*100);
end

figure;
prior_precis_colors = jet(2*length(pp_bins)-1);
prior_precis_colors = flipud(prior_precis_colors(1:length(pp_bins)-1,:));
for i = 1:size(prior_precis_colors,1)
    lineProps.col{i} = prior_precis_colors(i,:);
end
mseb(sp_axis,latency_means',latency_errors',lineProps)

xlim([sp_axis(1) sp_axis(end)])
xlabel('Sensory uncertainty','FontSize',20)
ylabel('Average decision latency (timesteps)','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
title('Prior Precision vs. Sensory Uncertainty: Decision Latency','FontSize',24)


%% Now look at evolution of characteristic quantities (epistemic vs. instrumental value) within trials, see how it affects decision bound/latency

results_dir = uigetdir();

results_fnam = fullfile(results_dir,uigetfile('*.mat'));


%% load data 

load(results_fnam)

