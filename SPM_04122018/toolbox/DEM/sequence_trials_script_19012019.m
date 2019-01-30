
% Use a super high chi-value, combined with really low precisions (high sensory noise), to hopefully catch the 'linear part'
% of the increase of the posterior negentropy (stored as H_int in the MDP
% output structure), and approximate the drift diffusion style
% race-to-bound dynamics

results_dir = uigetdir();
addpath(genpath('/Users/conorheins/Documents/MATLAB/spm12_r6906_mexmaci64/'))

%% set up parameters
chi_val = 1.36; 

% noise values for parameterizing likelihood precision (these are really
% low compared to earlier simulations)
noise_values = exp(linspace(log(0.05),log(0.8),10)); 

prior_precision = 0.25; % confidence/precision of prior beliefs
N = 100; % number of trials per chi-value/precision-value pair

[ all_configs,scene_idx ] = generate_scenes();

prior_scene_belief = 'UP_RIGHT'; % just choose one prior_scene_belief for simplicity, can extend to full randomization later
true_scene = 14; % just choose one true_scene for simplicity, can extend to full randomization later (1 - 12 corresponds to UR, etc.)

T = 32; % number of temporal updates (sequential samples) at lower level
policy_depth = 32; % depth of policy evaluation for calculating expected free energy at lower level


results_array = [];

%% run simulations
for p_i = 1:length(noise_values)
    
    % randomly initialize agents' beliefs
%     prior_scene_beliefs = cell(N,1);
%     idx_tmp = randi([1 4],N,1);
%     prior_scene_beliefs(idx_tmp == 1) = {'UP_RIGHT'};
%     prior_scene_beliefs(idx_tmp == 2) = {'RIGHT_DOWN'};
%     prior_scene_beliefs(idx_tmp == 3) = {'DOWN_LEFT'};
%     prior_scene_beliefs(idx_tmp == 4) = {'LEFT_UP'};
%     clear idx_temp;

    
    % randomly initialize true hidden states
%     true_scenes = randi([1 48],N,1);
    
    for trial = 1:N
        
%         MDPdeep = initialize_MDPdeep(prior_scene_beliefs{trial},prior_precision,all_configs,scene_idx,true_scenes(trial));
        MDPdeep = initialize_MDPdeep(prior_scene_belief,prior_precision,all_configs,scene_idx,true_scene);
        MDPshallow = initialize_MDPshallow_v5(noise_values(p_i),T,policy_depth);
%         MDPshallow.chi = chi_values(chi_i);
        MDPshallow.chi = chi_val;
        % nest the shallow MDP within the deeper one
        MDPdeep.MDP = MDPshallow;
        MDPdeep.link = sparse(1,1,1,numel(MDPshallow.D),numel(MDPdeep.A));
        
        MDPfull(1,trial) = spm_MDP_check(MDPdeep);
        clear MDPdeep MDPshallow;
        
    end
    
    %% Solving
    
    % solve a sequence of trials
    %==========================================================================
    MDPresult  = spm_MDP_VB_X(MDPfull); clear MDPfull;
    
    % gather trial-wise info (time-course of H_int as well
                        
    for trial = 1:N
        
        curr_trial = MDPresult(trial);
        
        quadrant_saccades = find(and(curr_trial.o(1,:) > 1,curr_trial.o(1,:) < 5));
        
        trial_array = [];
        
        for ii = 1:length(quadrant_saccades)
            
            sacc_ii = curr_trial.mdp(quadrant_saccades(ii));
            
            sub_array = zeros(2,T + 1); % one extra column (on the left) to store saccade number
            
            sub_array(:,1) = quadrant_saccades(ii);
            sub_array(1,2:(length(sacc_ii.H_int)+1)) = sacc_ii.H_int;
            sub_array(2,2:(length(sacc_ii.rt)+1)) = sacc_ii.rt;
            trial_array = [trial_array; [repmat(trial,2,1), sub_array]];
            
        end
        
        results_array = [results_array;[repmat(noise_values(p_i),size(trial_array,1),1),trial_array]];
        
    end

    save(fullfile(results_dir,sprintf('RDPsearch_Hierarch11_ChiId%d_PId%d.mat',1,p_i)),'MDPresult','chi_val','noise_values','T','policy_depth',...
        'all_configs','scene_idx','prior_scene_belief','true_scene','-v7.3')

%     save(fullfile(results_dir,sprintf('RDPsearch_Hierarch11_ChiId%d_PId%d.mat',1,p_i)),'MDPresult','chi_val','noise_values','T','policy_depth',...
%         'all_configs','scene_idx','prior_scene_belief','true_scene','-v7.3')
    
    clear MDPresult curr_trial quadrant_saccades trial_array sacc_ii sub_array trial ii
    
end

save(fullfile(results_dir,'RDPsearch_Hierarch11_ChiId1_ResultsArray.mat'),'chi_val','results_array','noise_values','-v7.3');

clear all; close all;
%% accumulate results_array saved data, if results array hasn't already been accumulated within simulation loop

results_dir = uigetdir();

T = 32;

data_files = dir([results_dir,filesep,'*.mat']);
data_files = {data_files(:).name};

results_array = [];

for f_i = 1:length(data_files)
    load(fullfile(results_dir,data_files{f_i}));
    
    N = length(MDPresult);
    
    noise_idx = str2double(data_files{f_i}(strfind(data_files{f_i},'PId')+3));
    
    noise_val = noise_values(noise_idx);
    
    for trial = 1:N
        
        curr_trial = MDPresult(trial);
        
        quadrant_saccades = find(and(curr_trial.o(1,:) > 1,curr_trial.o(1,:) < 5));
        
        trial_array = [];
        
        for ii = 1:length(quadrant_saccades)
            
            sacc_ii = curr_trial.mdp(quadrant_saccades(ii));
            
            sub_array = zeros(2,T + 1); % one extra column (on the left) to store saccade number
            
            sub_array(:,1) = quadrant_saccades(ii);
            sub_array(1,2:(length(sacc_ii.H_int)+1)) = sacc_ii.H_int;
            sub_array(2,2:(length(sacc_ii.rt)+1)) = sacc_ii.rt;
            trial_array = [trial_array; [repmat(trial,2,1), sub_array]];
            
        end
        
        results_array = [results_array;[repmat(noise_val,size(trial_array,1),1),trial_array]];
    end
end


%% ~analysis mode~

results_dir = uigetdir();
load(fullfile(results_dir,uigetfile('*.mat')))

all_negEnt = results_array(1:2:end,:);
all_rt = results_array(2:2:end,:);

breakearly_idx = [];
time2break = [];

for ii = 1:size(all_negEnt,1)
    
    if any(all_negEnt(ii,4:end)==0)
        breakearly_idx = [breakearly_idx,ii];
        time2break = [time2break; find(all_negEnt(ii,4:end)==0,1)];
    end
    
end


breakearly_negEnt = all_negEnt(breakearly_idx,:);
breakearly_rt = all_rt(breakearly_idx,:);

latebreak_negEnt = breakearly_negEnt(time2break > 2,:);
latebreak_rt = breakearly_rt(time2break  > 2,:);

noise_colors = cool(length(noise_values));

for ii = 1:length(noise_values)
  
    
%     data2plot = latebreak_negEnt(latebreak_negEnt(:,1)==noise_values(ii),4:end).*latebreak_rt(latebreak_negEnt(:,1)==noise_values(ii),4:end);
% %     plot(data2plot','Color',noise_colors(ii,:));
%     plot(mean(data2plot,1),'DisplayName',sprintf('Sensory precision: %.2f',noise_values(ii)));
      
    plot(mean(latebreak_negEnt(latebreak_negEnt(:,1)==noise_values(ii),4:end),1),...
        'DisplayName',sprintf('Sensory precision: %.2f',noise_values(ii)))
    
    
%     plot(latebreak_negEnt(latebreak_negEnt(:,1)==noise_values(ii),4:end)','Color',noise_colors(ii,:));
    
    pause;
    
    hold on;
    
end
legend('show')



reaction_times = time2break .* sum(all_rt(breakearly_idx,4:end),2); % scale time taken to reach Occam's bound by the reaction time incurred
scatter(all_rt(breakearly_idx,1),reaction_times);

negEnt_woEB = all_negEnt;
negEnt_woEB(breakearly_idx,:) = [];

rt_woEB = all_rt;
rt_woEB(breakearly_idx,:) = [];

figure;

noise_colors = cool(length(noise_values));

increasing_flag = negEnt_woEB(:,10) > negEnt_woEB(:,4); % flag for neg-entropy increasing over time (decreasing entropy)
decreasing_flag = ~increasing_flag; % flag for decreasing neg-entropy over time (increasing entropy?)
beginning_of_trial_flag = negEnt_woEB(:,3) <= 5;

steepness = [];

for ii = 1:length(noise_values)
    
   negEnt_this_ii = negEnt_woEB(and(and(negEnt_woEB(:,1)==noise_values(ii),increasing_flag),beginning_of_trial_flag),:);    
%    plot(negEnt_this_ii(:,4:end)','Color',noise_colors(ii,:));

    mean_trajectory = mean(negEnt_this_ii(:,4:end),1);
    mean_trajectory = mean_trajectory - mean_trajectory(1); % start them all at same vertical shift
    plot(mean_trajectory,'Color',noise_colors(ii,:),'DisplayName',sprintf('Sensory precision: %.2f',noise_values(ii)));
%    plot(-abs(diff(negEnt_this_ii(:,4:end),1,2))','Color',noise_colors(ii,:));
%    abs_differences = mean(abs(diff(negEnt_this_ii(:,4:end),1,2)),1);
%    plot(abs_differences,'Color',noise_colors(ii,:));

   hold on;
   
   pause; 
   
%    steepness = [steepness; [repmat(noise_values(ii),size(negEnt_this_ii,1),1),negEnt_this_ii(:,10)-negEnt_this_ii(:,5)]];
    
end

legend('show')



    
    
    
    
    

