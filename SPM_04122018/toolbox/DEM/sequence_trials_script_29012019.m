% Use a super high chi-value, combined with really low precisions (high sensory noise), to hopefully catch the 'linear part'
% of the increase of the posterior negentropy (stored as H_int in the MDP
% output structure), and approximate the drift diffusion style
% race-to-bound dynamics

results_dir = uigetdir();
addpath(genpath('/Users/conorheins/Documents/MATLAB/spm12_r6906_mexmaci64/'))

%% set up parameters
chi_val = 1/10e30; 

% noise values for parameterizing likelihood precision (these are really
% low compared to earlier simulations)
noise_values = exp(linspace(log(0.1),log(3),25)); 

prior_precision = 0.25; % confidence/precision of prior beliefs
N = 50; % number of trials per chi-value/precision-value pair

[ all_configs,scene_idx ] = generate_scenes();

prior_scene_belief = 'UP_RIGHT'; % just choose one prior_scene_belief for simplicity, can extend to full randomization later
true_scene = 14; % just choose one true_scene for simplicity, can extend to full randomization later (1 - 12 corresponds to UR, etc.)

T = 20; % number of temporal updates (sequential samples) at lower level
policy_depth = 1; % depth of policy evaluation for calculating expected free energy at lower level

results_array = [];

%% run simulations
for p_i = 1:length(noise_values)

    
    for trial = 1:N
        
        MDPdeep = initialize_MDPdeep(prior_scene_belief,prior_precision,all_configs,scene_idx,true_scene);
        MDPshallow = initialize_MDPshallow_v7(noise_values(p_i),T,policy_depth);
        MDPshallow.chi = chi_val;
        MDPdeep.MDP = MDPshallow;
        MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));
        
        MDPfull(1,trial) = spm_MDP_check(MDPdeep);
        clear MDPdeep MDPshallow;
        
    end
    
    %% Solving
    
    % solve a sequence of trials
    %==========================================================================
    MDPresult  = spm_MDP_VB_X_RCH(MDPfull); clear MDPfull;
    
    % gather trial-wise info (time-course of H_int as well
                        
    for trial = 1:N
        
        curr_trial = MDPresult(trial);
        
        quadrant_saccades = find(and(curr_trial.o(1,:) > 1,curr_trial.o(1,:) < 6));
        
        trial_array = [];
        
        for ii = 1:length(quadrant_saccades)
            
            sacc_ii = curr_trial.mdp(quadrant_saccades(ii));
            
            sub_array = zeros(1,2); % one extra column (on the left) to store saccade number
            
            sub_array(1,1) = quadrant_saccades(ii);
            
            break_time = size(sacc_ii.o,2); %the first time they made a decision at the lower level ('perceptual policy')
            sub_array(1,2) = break_time;

            trial_array = [trial_array; [trial, sub_array]];
            
        end
        
        results_array = [results_array;[repmat(noise_values(p_i),size(trial_array,1),1),trial_array]];
        
    end

    save(fullfile(results_dir,sprintf('RDPsearch_Hierarch14_ChiId%d_PId%d.mat',1,p_i)),'MDPresult','chi_val','noise_values','T','policy_depth',...
        'all_configs','scene_idx','prior_scene_belief','true_scene','-v7.3')
    
    clear MDPresult curr_trial quadrant_saccades trial_array sacc_ii break_time sub_array trial ii
    
end

save(fullfile(results_dir,'RDPsearch_Hierarch14_ChiId1_ResultsArray.mat'),'chi_val','results_array','noise_values','-v7.3');

clear all; close all;


%% ~analysis mode~

addpath(genpath('/Users/conorheins/Documents/MATLAB/spm12_r6906_mexmaci64/'))

results_dir = uigetdir();
load(fullfile(results_dir,uigetfile('*.mat')))

stats_table = zeros(2,length(noise_values));
for p_i = 1:length(noise_values)
    
    stats_table(1,p_i) = mean(results_array(results_array(:,1)==noise_values(p_i),4));
    stats_table(2,p_i) = std(results_array(results_array(:,1)==noise_values(p_i),4))./sqrt(sum(results_array(:,1)==noise_values(p_i)));
    
end

figure(1)
barwitherr(stats_table(2,:),stats_table(1,:));
xticks([1 5 10 15 20 25]);
xtick_label_formatter = @(flt)sprintf('%.1f',flt);
xtick_names = cellfun(@(flt) xtick_label_formatter(flt),num2cell(noise_values([1 5 10 15 20 25]))','UniformOutput', false);
xticklabels(xtick_names);
xlabel('Sensory precision');
ylabel('Latency');

%% create evidence accumulation curves (classic ones, in terms of log Bayes
% factors) 

evidence4truestates = [];
evidence4otherstates = [];


for p_i = 1:length(noise_values)
    
    load(fullfile(results_dir,sprintf('RDPsearch_Hierarch14_ChiId1_PId%d.mat',p_i)));
    
    sub_results = results_array(results_array(:,1) == noise_values(p_i),:);
    
    for trial = 1:length(MDPresult)
        
        trial_results = sub_results(sub_results(:,2) == trial,:);
        
        sacc_idx = trial_results(:,3);
        
        for sacc_i = 1:length(sacc_idx)
                
            % create log Bayes factor curves for each alternative and plot
            % decision boundaries, given when the agent broke
            
            time2break = trial_results(trial_results(:,3)==sacc_idx(sacc_i),4);
            
            obs = MDPresult(trial).mdp(sacc_idx(sacc_i)).o(1,1:time2break-1) - 1; 
            
            A = MDPresult(trial).mdp(sacc_idx(sacc_i)).A{1}(2:end,2:end,1);
            num_outcomes = size(A,1);
            
            evidence_vecs = -ones(num_outcomes)/(num_outcomes-1) + ((num_outcomes)/(num_outcomes-1))*eye(num_outcomes);

            x = zeros(num_outcomes,time2break-1);
            
            for s_i = 1:num_outcomes
                x(s_i,1) = evidence_vecs(s_i,:)*log(A(obs(1),:)');
            end
            
            for s_i = 1:num_outcomes
                for t = 2:time2break-1
                    x(s_i,t) = x(s_i,t-1) + evidence_vecs(s_i,:) * log(A(obs(t),:)');
                end
            end
            
            true_state_idx = MDPresult(trial).mdp(sacc_idx(sacc_i)).s(1,1)-1;
            ts_vec = zeros(1,T+2);
            ts_vec(1) = noise_values(p_i);
            ts_vec(2) = sacc_idx(sacc_i);
            ts_vec(3:(3+time2break-2)) = x(true_state_idx,:);
            
            [~,os_idx] = setdiff(1:num_outcomes,true_state_idx);
            os_vecs = zeros(length(os_idx),T+2);
            os_vecs(:,1) = noise_values(p_i);
            os_vecs(:,2) = sacc_idx(sacc_i);
            os_vecs(:,3:(3+time2break-2)) = x(os_idx,:);
            
            evidence4truestates = [evidence4truestates;ts_vec];
            evidence4otherstates = [evidence4otherstates;os_vecs];            
            
        end
        
    end
    
end

%% do some plotting

break_times_ts = zeros(size(evidence4truestates,1),1);
break_times_os = zeros(size(evidence4otherstates,1),1);

for i = 1:size(evidence4truestates,1)
    break_times_ts(i) = find(evidence4truestates(i,:)==0,1)-2;
end

for i = 1:size(evidence4otherstates,1)
    break_times_os(i) = find(evidence4otherstates(i,:)==0,1)-2;
end

break_times_unique = unique(break_times_ts);


for i = 2:length(break_times_unique)
    
    temp_all = evidence4truestates(break_times_ts==break_times_unique(i),:);
    
    unique_p = unique(temp_all(:,1));
    
    means_std_temp = zeros(length(unique_p),T,2);
    for j = 1:length(unique_p)
        trajs = temp_all(temp_all(:,1)==unique_p(j),3:end);
        means_std_temp(j,:,1) = mean(trajs,1);
        means_std_temp(j,:,2) = std(trajs,0,1)./sqrt(size(trajs,1));
    end
    
    uncert_colors = jet(ceil(1.5*length(unique_p)));
    uncert_colors = flipud(uncert_colors(1:length(unique_p),:));
    for c_i = 1:size(uncert_colors,1)
        lineProps.col{c_i} = uncert_colors(c_i,:);
    end
    
    mseb(1:break_times_unique(i)-1,means_std_temp(:,1:break_times_unique(i)-1,1),means_std_temp(:,1:break_times_unique(i)-1,2),...
        lineProps);
    
    uncertain_labels = cell(1,length(unique_p));
    for c_i = 1:length(uncertain_labels)
        uncertain_labels{1,c_i} = sprintf('%.2f sensory precision',unique_p(c_i));
    end
    
    xlim([1 break_times_unique(i)-1])
    xlabel('Sampling index','FontSize',20)
    ylabel('log Evidence','FontSize',20)
    xt = get(gca,'XTickLabel');
    yt = get(gca,'YTickLabel');
    set(gca,'XTickLabel',xt,'fontsize',16);
    set(gca,'YTickLabel',yt,'fontsize',16);
    legend(uncertain_labels,'Location','bestoutside','FontSize',16)
    
    pause; close gcf;
    
end


%% latency distributions

for p_i = 16:length(noise_values)
    
    [n,edges] = histcounts(results_array(results_array(:,1)==noise_values(p_i),4),10,'Normalization','probability');
    
    x_axis = edges(1:end-1) + diff(edges)./2;
    hold on; plot(x_axis,smooth(n,5),'LineWidth',1,'DisplayName',sprintf('Sensory Precision: %.1f',noise_values(p_i)));
    
end
legend('show')




