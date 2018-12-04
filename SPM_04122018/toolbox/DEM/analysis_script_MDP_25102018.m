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
load('stats_matrix_251018.mat')
disp(var_names)

%% Compare trajectories of epistemic value for different uncertainty and prior precision values

% look at trajectories, only for trials with constant number of timesteps
% and for one prior precision

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));

% u_bins = linspace(min(uncertainties),max(uncertainties),11);

u_bins = linspace(log(min(uncertainties)),log(max(uncertainties)),5);
u_bins = logspace(min(uncertainties),max(uncertainties),5);
pp_bins = linspace(min(prior_precisions),max(prior_precisions),4);

trial_length = 5;
trial_length_idx = master_stats_matrix(:,strcmp(var_names,'Time to Choose')) == trial_length;

EV_columns = find(strcmp(var_names,'EV1')):find(strcmp(var_names,sprintf('EV%d',trial_length)));

for pp_i = 1:length(pp_bins)-1
    
    p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(pp_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(pp_i+1));
    this_pbin_data = zeros(length(u_bins)-1,2,length(EV_columns));
    for u_i = 1:length(u_bins)-1
        u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1));
        trajs = master_stats_matrix(trial_length_idx&p_bin_idx&u_bin_idx,EV_columns);
        this_pbin_data(u_i,1,:) = mean(trajs,1);
        this_pbin_data(u_i,2,:) = std(trajs,0,1)./sqrt(size(trajs,1));
    end
    
    means = squeeze(this_pbin_data(:,1,:));
    errors = squeeze(this_pbin_data(:,2,:));
    
    
    figure()
    uncert_colors = jet(ceil(1.5*length(u_bins)-1));
    uncert_colors = flipud(uncert_colors(1:length(u_bins)-1,:));
    for i = 1:size(uncert_colors,1)
        lineProps.col{i} = uncert_colors(i,:);
    end

    mseb(1:trial_length,means,errors,lineProps)

    uncertain_labels = cell(1,length(u_bins)-1);
    for i = 1:length(uncertain_labels)
        uncertain_labels{1,i} = sprintf('%.2f - %.2f sensory noise level',u_bins(i),u_bins(i+1));
    end
    
    xlim([1 trial_length])
    xlabel('Time before categorization','FontSize',20)
    ylabel('Epistemic Value','FontSize',20)
    xt = get(gca,'XTickLabel');
    yt = get(gca,'YTickLabel');
    set(gca,'XTickLabel',xt,'fontsize',16);
    set(gca,'YTickLabel',yt,'fontsize',16);
    legend(uncertain_labels,'Location','bestoutside','FontSize',16)
    title(sprintf('Epistemic Value over Time for different Sensory Uncertainties: Prior Precision %.0f - %.0f %%',100*pp_bins(pp_i),100*pp_bins(pp_i+1)),'FontSize',24)
    
    saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/',sprintf('EVTrajectory_triallength%d_HighpriorPrecis.png',trial_length)));

end

%% now do the same, but binning across 3 different uncertainty levels 
% (one figure per uncertainty level), and showing different priors on the same
% plot

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));

u_bins = linspace(min(uncertainties),max(uncertainties),4);
pp_bins = linspace(min(prior_precisions),max(prior_precisions),5);

trial_length = 4;
trial_length_idx = master_stats_matrix(:,strcmp(var_names,'Time to Choose')) == trial_length;

EV_columns = find(strcmp(var_names,'EV1')):find(strcmp(var_names,sprintf('EV%d',trial_length)));

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1));
    this_ubin_data = zeros(length(pp_bins)-1,2,length(EV_columns));
    for p_i = 1:length(pp_bins)-1
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(p_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(p_i+1));
        trajs = master_stats_matrix(trial_length_idx&u_bin_idx&p_bin_idx,EV_columns);
        this_ubin_data(p_i,1,:) = mean(trajs,1);
        this_ubin_data(p_i,2,:) = std(trajs,0,1)./sqrt(size(trajs,1));
    end
    
    means = squeeze(this_ubin_data(:,1,:));
    errors = squeeze(this_ubin_data(:,2,:));
    
    figure;
    prior_precis_colors = jet(2*length(pp_bins)-1);
    prior_precis_colors = flipud(prior_precis_colors(1:length(pp_bins)-1,:));
    for i = 1:size(prior_precis_colors,1)
        lineProps.col{i} = prior_precis_colors(i,:);
    end

    mseb(1:trial_length,means,errors,lineProps)

    prior_precision_labels = cell(1,length(pp_bins)-1);
    for i = 1:length(prior_precision_labels)
        prior_precision_labels{1,i} = sprintf('%.0f - %0.f %% prior belief',pp_bins(i)*100,pp_bins(i+1)*100);
    end
    
    xlim([1 trial_length])
    xlabel('Time before categorization','FontSize',20)
    ylabel('Epistemic Value','FontSize',20)
    xt = get(gca,'XTickLabel');
    yt = get(gca,'YTickLabel');
    set(gca,'XTickLabel',xt,'fontsize',16);
    set(gca,'YTickLabel',yt,'fontsize',16);
    legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
    title(sprintf('Epistemic Value over Time for different Prior Precisions: Sensory Uncertainty %.2f - %.2f %%',u_bins(u_i),u_bins(u_i+1)),'FontSize',24)
    
    saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/',sprintf('EVTrajectory_triallength%d_HighNoise.png',trial_length)));

end



%% averaged only for trials with the same number of timesteps

trial_length = 5;

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
% u_bins = linspace(min(uncertainties),max(uncertainties),16);  
u_bins = linspace(min(uncertainties),max(uncertainties),6);
prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
pp_bins = linspace(min(prior_precisions),max(prior_precisions),16); 
EV_data = zeros(length(u_bins)-1,length(pp_bins)-1,2);

trial_length_idx = master_stats_matrix(:,strcmp(var_names,'Time to Choose')) == trial_length;

EV_columns = find(strcmp(var_names,'EV1')):find(strcmp(var_names,sprintf('EV%d',trial_length)));

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1)); 
    for pp_i = 1:length(pp_bins)-1
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(pp_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(pp_i+1));
        EV_average = mean(master_stats_matrix(trial_length_idx&u_bin_idx&p_bin_idx,EV_columns),2); % average of epistemic value over time, for timesteps before categorization
        EV_data(u_i,pp_i,1) = mean(EV_average); 
        EV_data(u_i,pp_i,2) = std(EV_average)./sqrt(length(EV_average));
    end
       
end

EV_means = squeeze(EV_data(:,:,1));
EV_errors = squeeze(EV_data(:,:,2));
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
mseb(uncertainty_axis,EV_means',EV_errors',lineProps)

xlim([uncertainty_axis(1) uncertainty_axis(end)])
xlabel('Sensory uncertainty','FontSize',20)
ylabel('Time-averaged Epistemic Value','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
title(sprintf('Prior Precision vs. Sensory Uncertainty: Average Epistemic Value for %d-saccade Trials',trial_length),'FontSize',24)

% saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/',sprintf('AverageEpValue_triallength%d.png',trial_length)));


%% make scatter plots of sensory uncertainty vs. epistemic value, for different prior precisions

% separate_by = 'belief_coherence';
separate_by = [];
trial_length = 2;

low_precis_top = 0.4; % highest bin-edge of the 'low precision' agents
mid_precis_top = 0.7; % highest bin-edge of the 'medium precision' agents

labels{1} = sprintf('Prior Precision: 25 %% to %0.f %%',low_precis_top*100);
labels{2} = sprintf('Prior Precision: %0.f %% to %0.f %%',(low_precis_top*100)+1,mid_precis_top*100);

if isempty(separate_by)
    figure(1);
    
    if mid_precis_top <= max(prior_precisions)
        labels{3} = sprintf('Prior Precision: %0.f %% to %0.f %%',(mid_precis_top*100)+1,max(prior_precisions)*100);
    end
    EV_columns = find(strcmp(var_names,'EV1')):find(strcmp(var_names,sprintf('EV%d',trial_length)));
    
    trial_length_idx = master_stats_matrix(:,strcmp(var_names,'Time to Choose')) == trial_length;
    
    uncertainties = master_stats_matrix(trial_length_idx,strcmp(var_names,'Uncertainty'));
    prior_precisions = master_stats_matrix(trial_length_idx,strcmp(var_names,'Prior Belief Precision'));
    averageEV = mean(master_stats_matrix(trial_length_idx,EV_columns),2);
    
    
    % low precision scatter plots
    lowPrecis = prior_precisions <= low_precis_top;
    colors_lowPrecis = repmat([0 0 .8],length(find(lowPrecis)),1);
   
    scatter(uncertainties(lowPrecis),averageEV(lowPrecis),30,colors_lowPrecis,'filled');
    hold on;
    
    % mid precision scatter plots
    midPrecis = and(prior_precisions > low_precis_top,prior_precisions <= mid_precis_top);
    colors_midPrecis = repmat([0 0.7 .5],length(find(midPrecis)),1);
   
    scatter(uncertainties(midPrecis),averageEV(midPrecis),30,colors_midPrecis,'filled')
   
    if length(labels) == 3
        % high precision scatter plots
        highPrecis = prior_precisions > mid_precis_top;
        colors_highPrecis = repmat([0.6 0.1 0.4],length(find(highPrecis)),1);
        
        scatter(uncertainties(highPrecis),averageEV(highPrecis),30,colors_highPrecis,'filled');
    end

    legend(labels,'Location','bestoutside','FontSize',16);
    
    xlim([min(uncertainties) max(uncertainties)])
    xlabel('Sensory uncertainty','FontSize',20)
    ylabel(sprintf('Epistemic Value averaged over %d timesteps',trial_length),'FontSize',20)
    xt = get(gca,'XTickLabel');
    yt = get(gca,'YTickLabel');
    set(gca,'XTickLabel',xt,'fontsize',16);
    set(gca,'YTickLabel',yt,'fontsize',16);
    title(sprintf('Prior Precision vs. Sensory Uncertainty: Average Epistemic Value for %d-saccade Trials',trial_length),'FontSize',24)
    
    clear labels;


elseif strcmp(separate_by,'belief_coherence')
    
    figure(1);
    
    labels_correct_only = {'Low Precision, Correct Beliefs','Mid Precision, Correct Beliefs','High Precision, Correct Beliefs'};
    labels_all = {'Low Precision, Correct Beliefs','Mid Precision, Correct Beliefs','High Precision, Correct Beliefs','Low Precision: Incorrect Beliefs','Mid Precision, Incorrect Beliefs',...
        'High Precision, Incorrect Beliefs'};
    
    EV_columns = find(strcmp(var_names,'EV1')):find(strcmp(var_names,sprintf('EV%d',trial_length)));
    
    trial_length_idx = master_stats_matrix(:,strcmp(var_names,'Time to Choose')) == trial_length;
    
    true_scenes = master_stats_matrix(trial_length_idx,strcmp(var_names,'True Scene'));
    scene_beliefs = master_stats_matrix(trial_length_idx,strcmp(var_names,'Scene Belief'));
    uncertainties = master_stats_matrix(trial_length_idx,strcmp(var_names,'Uncertainty'));
    prior_precisions = master_stats_matrix(trial_length_idx,strcmp(var_names,'Prior Belief Precision'));
    averageEV = mean(master_stats_matrix(trial_length_idx,EV_columns),2);
    
   
    % low precision scatter plots
    lowPrecis_correct = prior_precisions <= low_precis_top & true_scenes == scene_beliefs;
    colors_lowPrecis_correct = repmat([0 0 .8],length(find(lowPrecis_correct)),1);
    lowPrecis_incorrect = prior_precisions <= low_precis_top & true_scenes ~= scene_beliefs;
    colors_lowPrecis_incorrect = repmat([0 0.6 1],length(find(lowPrecis_incorrect)),1);
    
    scatter(uncertainties(lowPrecis_correct),averageEV(lowPrecis_correct),30,colors_lowPrecis_correct,'filled');
    hold on;
    
    % mid precision scatter plots
    midPrecis_correct = and(prior_precisions > low_precis_top,prior_precisions <= mid_precis_top) & true_scenes == scene_beliefs;
    colors_midPrecis_correct = repmat([0 0.7 .5],length(find(midPrecis_correct)),1);
    midPrecis_incorrect = and(prior_precisions > low_precis_top,prior_precisions <= mid_precis_top) & true_scenes ~= scene_beliefs;
    colors_midPrecis_incorrect = repmat([.2 1 0.5],length(find(midPrecis_incorrect)),1);
    
    scatter(uncertainties(midPrecis_correct),averageEV(midPrecis_correct),30,colors_midPrecis_correct,'filled')
    
    % high precision scatter plots
    highPrecis_correct = prior_precisions > mid_precis_top & true_scenes == scene_beliefs;
    colors_highPrecis_correct = repmat([0.6 0.1 0.4],length(find(highPrecis_correct)),1);
    highPrecis_incorrect = prior_precisions > mid_precis_top & true_scenes ~= scene_beliefs;
    colors_highPrecis_incorrect = repmat([0.9 0.3 0.2],length(find(highPrecis_incorrect)),1);
    
    scatter(uncertainties(highPrecis_correct),averageEV(highPrecis_correct),30,colors_highPrecis_correct,'filled');
    
    legend(labels_correct_only,'Location','bestoutside','FontSize',16);

    %% now add the incorrect ones (across all precision levels)
    
    pause;
    
    scatter(uncertainties(lowPrecis_incorrect),averageEV(lowPrecis_incorrect),30,colors_lowPrecis_incorrect,'filled');
    
    scatter(uncertainties(midPrecis_incorrect),averageEV(midPrecis_incorrect),30,colors_midPrecis_incorrect,'filled')
    
    scatter(uncertainties(highPrecis_incorrect),averageEV(highPrecis_incorrect),30,colors_highPrecis_incorrect,'filled');
    
    
    legend(labels_all,'Location','bestoutside','FontSize',16);
end

%% save figure after adjusting
saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/',sprintf('AvgEV_%dsaccades.png',trial_length)));

%% number of quadrant revisit data, incorrect trials (with new format, no more revisits after categorization)

% agreement_idx = master_stats_matrix(:,strcmp(var_names,'True Scene')) == master_stats_matrix(:,strcmp(var_names,'Scene Belief'));
correct_idx = master_stats_matrix(:,strcmp(var_names,'Accuracy')) == 1;
incorrect_idx = master_stats_matrix(:,strcmp(var_names,'Accuracy')) == 2;
non_nan_idx = ~isnan(master_stats_matrix(:,strcmp(var_names,'Number Quadrant Revisits')));

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = linspace(min(uncertainties),max(uncertainties),16);  

prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
pp_bins = linspace(min(prior_precisions),max(prior_precisions),16); 

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

%%

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
        num_revisits = master_stats_matrix(~agreement_idx&u_bin_idx&p_bin_idx,strcmp(var_names,'Number Quadrant Revisits'));
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
title('Prior Precision vs. Sensory Uncertainty: Quadrant Revisits','FontSize',24)

%%

saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/','numquad_revisit_incorrect.png'))


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

%%
saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/','ambivalence.png'))


%% now do straight average epistemic value

agreement_idx = master_stats_matrix(:,strcmp(var_names,'True Scene')) == master_stats_matrix(:,strcmp(var_names,'Scene Belief'));

uncertainties = master_stats_matrix(:,strcmp(var_names,'Uncertainty'));
u_bins = linspace(min(uncertainties),max(uncertainties),16);  % better for looking at the interaction with incorrect prior beliefs

prior_precisions = master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision'));
pp_bins = linspace(min(prior_precisions),max(prior_precisions),16); % better for looking at the interaction with incorrect prior beliefs

EV_data = zeros(length(u_bins)-1,length(pp_bins)-1,2);

% trial_length = 7;
% trial_length_idx = master_stats_matrix(:,strcmp(var_names,'Time to Choose')) == trial_length;
% 
% EV_columns = find(strcmp(var_names,'EV1')):find(strcmp(var_names,sprintf('EV%d',trial_length)));
EV_columns = find(strcmp(var_names,'EV1')):find(strcmp(var_names,'EV8'));

for u_i = 1:length(u_bins)-1
    
    u_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Uncertainty')) >= u_bins(u_i),master_stats_matrix(:,strcmp(var_names,'Uncertainty')) < u_bins(u_i+1)); 
    for pp_i = 1:length(pp_bins)-1
        p_bin_idx = and(master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) >= pp_bins(pp_i),master_stats_matrix(:,strcmp(var_names,'Prior Belief Precision')) < pp_bins(pp_i+1));
        EV = mean(master_stats_matrix(agreement_idx&u_bin_idx&p_bin_idx,EV_columns),2);
        EV_data(u_i,pp_i,1) = mean(EV); 
        EV_data(u_i,pp_i,2) = std(EV)./sqrt(length(EV));
    end
       
end

EV_mean = squeeze(EV_data(:,:,1));
EV_errors = squeeze(EV_data(:,:,2));
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
mseb(uncertainty_axis,EV_mean',EV_errors',lineProps)

xlim([uncertainty_axis(1) uncertainty_axis(end)])
xlabel('Sensory uncertainty','FontSize',20)
ylabel('Time-averaged epistemic value','FontSize',20)
xt = get(gca,'XTickLabel');
yt = get(gca,'YTickLabel');
set(gca,'XTickLabel',xt,'fontsize',16);
set(gca,'YTickLabel',yt,'fontsize',16);
legend(prior_precision_labels,'Location','bestoutside','FontSize',16)
title('Prior Confidence vs. Sensory Uncertainty: Epistemic Value','FontSize',24)

%%
saveas(gca,fullfile('/Users/conorheins/Documents/Presentations/','EV_all_correct.png'))






    
