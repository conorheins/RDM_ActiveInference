%% analysis_script_MDP_171012018

% choose data with MDP structure to analyze

data_dir = 'C:/RDPsearch_Hierarch5_results/';

all_files = dir(fullfile(data_dir,'*.mat'));
all_files = {all_files(:).name};

num_files = length(all_files);

master_stats_matrix = [];

tic
for f_i = 1:num_files
    
    [sub_stats_matrix] = collect_stats_across_trialsv3(data_dir,all_files,f_i);
    
    master_stats_matrix = [master_stats_matrix;sub_stats_matrix];
    
    clear sub_stats_matrix;
    
end
fprintf('Time taken to load/analyze all %d files: %.2f minutes\n',num_files,toc/60)

var_names = {'Uncertainty', 'Prior Belief Precision', 'True Scene','Scene Belief','Realization Latency','Number Quadrant Revisits'...
    'Accuracy','dFE Value','dFE Time','Average LO','minLO','maxLO','EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','Time to Choose'}; 
save(fullfile(data_dir,sprintf('stats_matrix_%s.mat',datestr(now,'ddmmyy'))),'all_files','master_stats_matrix','var_names')

%% analyze stuff
load('masterstatsmatrix__new_25102018.mat')
disp(var_names)
% get the indices when the agent's belief cohered with the true belief
agreement_idx = master_stats_matrix(:,strcmp(var_names,'True Scene')) == master_stats_matrix(:,strcmp(var_names,'Scene Belief'));


%% time until 'eureka moment' 

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = linspace(min(uncertainties),max(uncertainties),16);  % better for looking at the interaction with incorrect prior beliefs
% u_bins = linspace(1.0,max(uncertainties),11); % better for looking at the interaction with correct prior beliefs

prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
pp_bins = linspace(min(prior_precisions),max(prior_precisions),16); % better for looking at the interaction with incorrect prior beliefs
% pp_bins = linspace(0.25,0.8,6); % better for looking at the interaction with correct prior beliefs

eureka_data = zeros(length(u_bins)-1,length(pp_bins)-1,2);

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1)); 
    for pp_i = 1:length(pp_bins)-1
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(pp_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(pp_i+1));
        time2realize = master_stats_matrix(~agreement_idx&u_bin_idx&p_bin_idx,strcmp(var_names,'Realization Latency'));
        eureka_data(u_i,pp_i,1) = mean(time2realize); 
        eureka_data(u_i,pp_i,2) = std(time2realize)./sqrt(length(time2realize));
    end
       
end

eureka_means = squeeze(eureka_data(:,:,1));
eureka_errors = squeeze(eureka_data(:,:,2));
uncertainty_axis = u_bins(1:end-1) + diff(u_bins)/2;
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
mseb(uncertainty_axis,eureka_means',eureka_errors',lineProps)

xlim([uncertainty_axis(1) uncertainty_axis(end)])
xlabel('Sensory uncertainty','FontSize',20)
ylabel('Timesteps until realization of true scene','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
title('Prior Precision vs. Sensory Uncertainty: Realization Latency','FontSize',24)

saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/','realization_latency.png'));

%% number of quadrant revisit data, incorrect trials (with new format, no more revisits after categorization)

% agreement_idx = master_stats_matrix(:,strcmp(var_names,'True Scene')) == master_stats_matrix(:,strcmp(var_names,'Scene Belief'));
correct_idx = master_stats_matrix(:,strcmp(var_names,'Accuracy')) == 1;
incorrect_idx = master_stats_matrix(:,strcmp(var_names,'Accuracy')) == 2;
non_nan_idx = ~isnan(master_stats_matrix(:,strcmp(var_names,'Number Quadrant Revisits')));

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = linspace(min(uncertainties),max(uncertainties),16);  % better for looking at the interaction with incorrect prior beliefs
% u_bins = linspace(min(uncertainties),2.78,16);
% u_bins = linspace(1.0,max(uncertainties),11); % better for looking at the interaction with correct prior beliefs

prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
pp_bins = linspace(min(prior_precisions),max(prior_precisions),16); % better for looking at the interaction with incorrect prior beliefs
% pp_bins = linspace(min(prior_precisions),0.77,16);
% pp_bins = linspace(0.25,0.8,6); % better for looking at the interaction with correct prior beliefs

revisit_data = zeros(length(u_bins)-1,length(pp_bins)-1,2);

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1)); 
    for pp_i = 1:length(pp_bins)-1
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(pp_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(pp_i+1));
        num_revisits = master_stats_matrix(non_nan_idx&incorrect_idx&u_bin_idx&p_bin_idx,strcmp(var_names,'Number Quadrant Revisits'));
        revisit_data(u_i,pp_i,1) = mean(num_revisits); 
        revisit_data(u_i,pp_i,2) = std(num_revisits)./sqrt(length(num_revisits));
    end
       
end

revisit_mean = squeeze(revisit_data(:,:,1));
revisit_errors = squeeze(revisit_data(:,:,2));
uncertainty_axis = u_bins(1:end-1) + diff(u_bins)/2;

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
mseb(uncertainty_axis,revisit_mean',revisit_errors',lineProps)

xlim([uncertainty_axis(1) uncertainty_axis(end)])
xlabel('Sensory uncertainty','FontSize',20)
ylabel('Number of quadrant revisits before choice','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
title('Prior Confidence vs. Sensory Uncertainty: epistemic foraging','FontSize',24)

saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/','numquad_revisit_incorrect.png'));

%% number of quadrant revisit data, correct trials (with new format, no more revisits after categorization)

agreement_idx = master_stats_matrix(:,strcmp(var_names,'True Scene')) == master_stats_matrix(:,strcmp(var_names,'Scene Belief'));

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = linspace(min(uncertainties),max(uncertainties),16);  % better for looking at the interaction with incorrect prior beliefs
% u_bins = linspace(1.0,max(uncertainties),11); % better for looking at the interaction with correct prior beliefs

prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
pp_bins = linspace(min(prior_precisions),max(prior_precisions),16); % better for looking at the interaction with incorrect prior beliefs
% pp_bins = linspace(0.25,0.8,6); % better for looking at the interaction with correct prior beliefs

revisit_data = zeros(length(u_bins)-1,length(pp_bins)-1,2);

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1)); 
    for pp_i = 1:length(pp_bins)-1
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(pp_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(pp_i+1));
        num_revisits = master_stats_matrix(agreement_idx&u_bin_idx&p_bin_idx,strcmp(var_names,'Number Quadrant Revisits'));
        revisit_data(u_i,pp_i,1) = mean(num_revisits); 
        revisit_data(u_i,pp_i,2) = std(num_revisits)./sqrt(length(num_revisits));
    end
       
end

revisit_mean = squeeze(revisit_data(:,:,1));
revisit_errors = squeeze(revisit_data(:,:,2));
uncertainty_axis = u_bins(1:end-1) + diff(u_bins)/2;

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
mseb(uncertainty_axis,revisit_mean',revisit_errors',lineProps)

xlim([uncertainty_axis(1) uncertainty_axis(end)])
xlabel('Sensory uncertainty','FontSize',20)
ylabel('Number of quadrant revisits before choice','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
title('Prior Confidence vs. Sensory Uncertainty: epistemic foraging','FontSize',24)

saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/','numquad_revisit_correct.png'))


%% Accuracy/Inaccuracy/Ambivalence as a function of prior beliefs and uncertainty

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = linspace(min(uncertainties),max(uncertainties),16);  % better for looking at the interaction with incorrect prior beliefs

prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
pp_bins = linspace(min(prior_precisions),max(prior_precisions),16); % better for looking at the interaction with incorrect prior beliefs

accur_data = zeros(length(u_bins)-1,length(pp_bins)-1,2);

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1)); 
    for pp_i = 1:length(pp_bins)-1
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(pp_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(pp_i+1));
        accuracy = master_stats_matrix(u_bin_idx&p_bin_idx,strcmp(var_names,'Accuracy')) == 3;
        accur_data(u_i,pp_i,1) = mean(accuracy); 
        accur_data(u_i,pp_i,2) = std(accuracy)./sqrt(length(accuracy));
    end
       
end

accur_mean = squeeze(accur_data(:,:,1));
accur_errors = squeeze(accur_data(:,:,2));
uncertainty_axis = u_bins(1:end-1) + diff(u_bins)/2;

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
mseb(uncertainty_axis,accur_mean',accur_errors',lineProps)

xlim([uncertainty_axis(1) uncertainty_axis(end)])
xlabel('Sensory uncertainty','FontSize',20)
ylabel('Proportion of uncategorized trials','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
title('Prior Confidence vs. Sensory Uncertainty: Ambivalence','FontSize',24)

saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/','ambivalence.png'))






    
