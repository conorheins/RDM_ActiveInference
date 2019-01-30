
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
noise_values = exp(linspace(log(0.4),log(0.8),11)); 

prior_precision = 0.25; % confidence/precision of prior beliefs
N = 50; % number of trials per chi-value/precision-value pair

[ all_configs,scene_idx ] = generate_scenes();

prior_scene_belief = 'UP_RIGHT'; % just choose one prior_scene_belief for simplicity, can extend to full randomization later
true_scene = 14; % just choose one true_scene for simplicity, can extend to full randomization later (1 - 12 corresponds to UR, etc.)

T = 32; % number of temporal updates (sequential samples) at lower level
policy_depth = 1; % depth of policy evaluation for calculating expected free energy at lower level

results_array = [];

%% run simulations
for p_i = 1:length(noise_values)

    
    for trial = 1:N
        
        MDPdeep = initialize_MDPdeep(prior_scene_belief,prior_precision,all_configs,scene_idx,true_scene);
        MDPshallow = initialize_MDPshallow_v6(noise_values(p_i),T,policy_depth);
        MDPshallow.chi = chi_val;
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
            
            sub_array = zeros(1,2); % one extra column (on the left) to store saccade number
            
            sub_array(1,1) = quadrant_saccades(ii);
            
            break_time = find(sacc_ii.o(3,:) > 1,1); %the first time they made a decision at the lower level ('perceptual policy')
            if ~isempty(break_time)
                sub_array(1,2) = break_time;
            else
                sub_array(1,2) = T;
            end
            trial_array = [trial_array; [trial, sub_array]];
            
        end
        
        results_array = [results_array;[repmat(noise_values(p_i),size(trial_array,1),1),trial_array]];
        
    end

    save(fullfile(results_dir,sprintf('RDPsearch_Hierarch12_ChiId%d_PId%d.mat',1,p_i)),'MDPresult','chi_val','noise_values','T','policy_depth',...
        'all_configs','scene_idx','prior_scene_belief','true_scene','-v7.3')
    
    clear MDPresult curr_trial quadrant_saccades trial_array sacc_ii break_time sub_array trial ii
    
end

save(fullfile(results_dir,'RDPsearch_Hierarch12_ChiId1_ResultsArray.mat'),'chi_val','results_array','noise_values','-v7.3');

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


for p_i = 1:length(noise_values)
    
    load(fullfile(results_dir,sprintf('RDPsearch_Hierarch12_ChiId1_PId%d.mat',p_i)));
    
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
            
            for s_i = 1:num_outcomes
                if MDPresult(trial).mdp(sacc_idx(sacc_i)).s(1,1)-1 == s_i
                    plot(x(s_i,:),'LineWidth',2);
                else
                    plot(x(s_i,:),'LineWidth',1);
                end
                
                hold on;
                
            end
            
            title(sprintf('Trial %d, Saccade No. %d, Precision %.2f',trial,sacc_idx(sacc_i),noise_values(p_i)));
            pause;            
            close gcf;
            
            
%             temp = squeeze(MDPresult(trial).mdp(sacc_idx(sacc_i)).xn{1}(end,2:end,:,:));
%             beliefs_squeezed = zeros(size(temp,1),size(temp,3));
%             for xn_i = 1:size(temp,1)
%                 beliefs_squeeze(xn_i,:) = diag(squeeze(temp(xn_i,:,:)));
%             end
%             
%             plot(beliefs_squeeze','LineWidth',2);
%             hold on; plot([trial_results(trial_results(:,3)==sacc_idx(sacc_i),4)-1, trial_results(trial_results(:,3)==sacc_idx(sacc_i),4)-1],...
%                 [0, 1],'r--','LineWidth',1);
%             
%             title(sprintf('Trial %d, Saccade No. %d, Precision %.2f',trial,sacc_idx(sacc_i),noise_values(p_i)));
%             pause;
%             
%             close gcf;
            
        end
        
    end
    
end
                
            
            
        
        
        
        
        
    
    
    

    

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



    
    
    
    
    

