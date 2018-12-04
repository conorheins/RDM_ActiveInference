%% analysis_script_MDP_211112018

% choose data with MDP structure to analyze

data_dir = 'C:/RDPsearch_Hierarch7_results/';

all_files = dir(fullfile(data_dir,'*.mat'));
all_files = {all_files(:).name};

num_files = length(all_files);

master_stats_matrix = [];

tic
for f_i = 1:num_files
    
    [sub_stats_matrix] = collect_stats_across_trialsv4(data_dir,all_files,f_i);
    
    master_stats_matrix = [master_stats_matrix;sub_stats_matrix];
    
    clear sub_stats_matrix;
    
end
fprintf('Time taken to load/analyze all %d files: %.2f minutes\n',num_files,toc/60)

var_names = {'Uncertainty', 'Prior Belief Precision', 'True Scene','Scene Belief','Realization Latency','Number Quadrant Revisits'...
    'Accuracy','dFE Value','dFE Time','Average LO','minLO','maxLO','Time to Choose','EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','Time to Bound'};
save(fullfile(data_dir,sprintf('stats_matrix_%s.mat',datestr(now,'ddmmyy'))),'all_files','master_stats_matrix','var_names')

%% analyze stuff
load stats_matrix_221118.mat

disp(var_names)

%% Compare race-to-bound times for different uncertainty and prior precision values

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = exp(linspace(log(min(uncertainties)),log(max(uncertainties)),25));

race2bound_data = zeros(length(u_bins)-1,2);

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1));
    temp = master_stats_matrix(u_bin_idx,strcmp(var_names,'Time to Bound'));
    race2bound_data(u_i,1) = mean(temp);
    race2bound_data(u_i,2) = std(temp)./sqrt(length(temp));
    %     race2bound_data(u_i,2) = std(temp);
    
end

% transformation uncertainty values into coherence values
uncertainty_axis = u_bins(1:end-1) + diff(u_bins)/2;
vectors = spm_softmax([uncertainty_axis;zeros(3,length(uncertainty_axis))]);
uncertainty_axis = vectors(1,:);

lineProps.col = {[0 0.5 0.25]};
figure;
mseb(uncertainty_axis,race2bound_data(:,1)',race2bound_data(:,2)',lineProps)

xlim([uncertainty_axis(1) uncertainty_axis(end)])
xlabel('% Coherence of fixated RDM','FontSize',20)
ylabel('Average time to reach bound for quadrant-visits','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
title('Race to Bound time as a function of coherence','FontSize',24)

%% save figure
saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/','AverageRace2Bound_vCoherence.png'));


%% plot trajectory of epistemic values, separated by prior precisions

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = exp(linspace(log(min(uncertainties)),log(max(uncertainties)),25));

prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
p_bins = exp(linspace(log(min(prior_precisions)),log(max(prior_precisions)),5));

EV_trajectories = zeros(length(u_bins)-1,length(p_bins)-1,8,2);

trial_length = 8;

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1));
    
    for p_i = 1:length(p_bins)-1
        
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= p_bins(p_i),...
            master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < p_bins(p_i+1));
        
        temp = master_stats_matrix(u_bin_idx&p_bin_idx,14:21);
        
        EV_trajectories(u_i,p_i,:,1) = mean(temp,1);
        EV_trajectories(u_i,p_i,:,2) = std(temp,0,1)./sqrt(size(temp,1));
        
    end
    
end

uncertainty_axis = u_bins(1:end-1) + diff(u_bins)/2;
vectors = spm_softmax([uncertainty_axis;zeros(3,length(uncertainty_axis))]);
uncertainty_axis = vectors(1,:);

for u_i = 1:length(uncertainty_axis)
    
    means = squeeze(EV_trajectories(u_i,:,:,1));
    errors = squeeze(EV_trajectories(u_i,:,:,2));
        
    figure;
    prior_precis_colors = jet(2*length(p_bins)-1);
    prior_precis_colors = flipud(prior_precis_colors(1:length(p_bins)-1,:));
    for i = 1:size(prior_precis_colors,1)
        lineProps.col{i} = prior_precis_colors(i,:);
    end

    mseb(1:trial_length,means,errors,lineProps)
    
    prior_precision_labels = cell(1,length(p_bins)-1);
    for i = 1:length(prior_precision_labels)
        prior_precision_labels{1,i} = sprintf('%.0f - %0.f %% prior belief',p_bins(i)*100,p_bins(i+1)*100);
    end
    
    xlim([1 trial_length])
    xlabel('Length of trial','FontSize',20)
    ylabel('Epistemic Value','FontSize',20)
    xt = get(gca,'XTickLabel');
    yt = get(gca,'YTickLabel');
    set(gca,'XTickLabel',xt,'fontsize',16);
    set(gca,'YTickLabel',yt,'fontsize',16);
    legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
    title(sprintf('Epistemic Value over Time for different Prior Precisions: RDM Coherence %.2f %%',uncertainty_axis(u_i) * 100),'FontSize',24)
    
    
    pause; close gcf;
end

%% Accuracy/Inaccuracy/Ambivalence as a function of prior beliefs and uncertainty

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = exp(linspace(log(min(uncertainties)),log(max(uncertainties)),6));
% u_bins = linspace(min(uncertainties),max(uncertainties),10);

prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
% p_bins = exp(linspace(log(min(prior_precisions)),log(max(prior_precisions)),11));
p_bins = [0.25 0.5 0.8];

accur_data = zeros(length(u_bins)-1,length(p_bins)-1,2);

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1)); 
    for pp_i = 1:length(p_bins)-1
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= p_bins(pp_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < p_bins(pp_i+1));
        accuracy = master_stats_matrix(u_bin_idx&p_bin_idx,strcmp(var_names,'Accuracy')) == 1;
        accur_data(u_i,pp_i,1) = mean(accuracy); 
        accur_data(u_i,pp_i,2) = std(accuracy)./sqrt(length(accuracy));
    end
       
end

accur_mean = squeeze(accur_data(:,:,1));
accur_errors = squeeze(accur_data(:,:,2));

uncertainty_axis = u_bins(1:end-1) + diff(u_bins)/2;
vectors = spm_softmax([uncertainty_axis;zeros(3,length(uncertainty_axis))]);
uncertainty_axis = vectors(1,:);

prior_precision_labels = cell(1,length(p_bins)-1);
for i = 1:length(prior_precision_labels)
    prior_precision_labels{1,i} = sprintf('%.0f - %0.f %% prior belief',p_bins(i)*100,p_bins(i+1)*100);
end

figure;
prior_precis_colors = jet(2*length(p_bins)-1);
prior_precis_colors = flipud(prior_precis_colors(1:length(p_bins)-1,:));
for i = 1:size(prior_precis_colors,1)
    lineProps.col{i} = prior_precis_colors(i,:);
end
mseb(uncertainty_axis,accur_mean',accur_errors',lineProps)

xlim([uncertainty_axis(1) uncertainty_axis(end)])
xlabel('RDM Coherence','FontSize',20)
ylabel('Proportion of uncategorized trials','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
title('Prior Confidence vs. Coherence: Accuracy','FontSize',24)

    
    
