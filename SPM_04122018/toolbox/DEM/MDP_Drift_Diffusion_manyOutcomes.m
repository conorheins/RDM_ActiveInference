function MDP_Drift_Diffusion_manyOutcomes(num_outcomes,precisions2test,T,num_repeats)

if nargin < 3 || ~exist('T','var') || isempty(T)
    T = 32; % default to 32 
end

if nargin < 4 || ~exist('num_repeats','var') || isempty(num_repeats)
    num_repeats = 1;
end

precis_colors = cool(length(precisions2test));

strings_true = @(string_i)(sprintf('Average trajectory for true state, precision: %.1f',string_i));
precis_labels_true = cellfun(@(x)strings_true(x),num2cell(precisions2test),'UniformOutput',false);

strings_other = @(string_i)(sprintf('Average trajectory for other states, precision: %.1f',string_i));
precis_labels_other = cellfun(@(x)strings_other(x),num2cell(precisions2test),'UniformOutput',false);

color_iter = 1;

all_means_true = zeros(length(precisions2test),T);
all_sems_true = zeros(length(precisions2test),T);

all_means_other = zeros(length(precisions2test),T);
all_sems_other = zeros(length(precisions2test),T);

for i = 1:length(precisions2test)
    precis = precisions2test(i);
%     D{1} = ones(num_outcomes,1);
    A{1} = spm_softmax(exp(precis)*log(eye(num_outcomes)+4));
%     B{1}(:,:,1) = eye(num_outcomes);           % Actions do nothing (the two B matrices here are only included to prevent the scheme entering 'HMM mode')
%     B{1}(:,:,2) = eye(num_outcomes);
%     C{1} = zeros(num_outcomes,1);
%     
%     for rep_i = 1:num_repeats
%         MDP(1,rep_i).A = A;
%         MDP(1,rep_i).B = B;
%         MDP(1,rep_i).C = C;
%         MDP(1,rep_i).D = D;
%         MDP(1,rep_i).T = T;
%         MDP(1,rep_i).s = randi(num_outcomes)*ones(1,T);
%     end
% 
%     mdp = spm_MDP_check(MDP);
%     MDP = spm_MDP_VB_X(mdp);


    true_states = randi(num_outcomes,num_repeats,1);
    
    x = zeros(num_outcomes,T,num_repeats);
%     all_beliefs = zeros(size(x));

    evidence_vecs = -ones(num_outcomes)/(num_outcomes-1) + ((num_outcomes)/(num_outcomes-1))*eye(num_outcomes);
    
    all_obs = zeros(1,T,num_repeats);
    for rep_i = 1:num_repeats
        
        obs = zeros(1,T);
        for j = 1:length(obs)
            obs(j) = find(rand < cumsum(A{1}(:,true_states(rep_i))),1);
        end
        
        for s_i = 1:num_outcomes
            x(s_i,1,rep_i) = evidence_vecs(s_i,:)*log(A{1}(obs(1),:)');
        end
        
        for s_i = 1:num_outcomes
            for t = 2:T
                x(s_i,t,rep_i) = x(s_i,t-1,rep_i) + evidence_vecs(s_i,:) * log(A{1}(obs(t),:)');
            end
        end
        
        all_obs(1,:,rep_i) = obs;
        
%         for s_i = 1:num_outcomes
%             x(s_i,1,rep_i) = evidence_vecs(s_i,:)*log(MDP(1,rep_i).A{1}(MDP(1,rep_i).o(1),:)');
%         end
%         
%         for s_i = 1:num_outcomes
%             for t = 2:T
%                 x(s_i,t,rep_i) = x(s_i,t-1,rep_i) + evidence_vecs(s_i,:) * log(MDP(1,rep_i).A{1}(MDP(1,rep_i).o(t),:)');
%             end
%         end
    
%         for s_i = 1:num_outcomes
%             all_beliefs(s_i,:,rep_i) = cumsum(diag(squeeze(MDP(1,rep_i).xn{1}(end,s_i,:,:))));
%         end
%         
%         all_beliefs(:,:,rep_i) = spm_softmax(squeeze(all_beliefs(:,:,rep_i)));
%         
%         negentropy = zeros(1,T,num_repeats);
%         for t = 1:MDP(1,rep_i).T
%             negentropy(1,t,rep_i) = all_beliefs(:,t,rep_i)'*log(all_beliefs(:,t,rep_i) + exp(-16));
%         end
%     
    end
    
    true_states_evidence = zeros(1,T,num_repeats);
    other_states_evidence = zeros(num_outcomes-1,T,num_repeats);
    
    for rep_i = 1:num_repeats
        
        state_idx = 1:num_outcomes;
%         [~,true_state_idx] = ismember(MDP(1,rep_i).s(1),state_idx);
        [~,true_state_idx] = ismember(true_states(rep_i),state_idx);
        
        true_states_evidence(1,:,rep_i) = x(true_state_idx,:,rep_i);
        other_states_evidence(:,:,rep_i) = x(state_idx(state_idx~=true_state_idx),:,rep_i);
        
    end
    
   
    if num_repeats <= 10
        for rep_i = 1:num_repeats
            plot(squeeze(true_states_evidence(1,:,rep_i)),'Color',precis_colors(color_iter,:),'LineWidth',0.1,'DisplayName',sprintf('Decision variable for true state, precision: %.1f',precis));
            hold on;
            plot(squeeze(sum(other_states_evidence(:,:,rep_i),1)),'Color',precis_colors(color_iter,:),'LineWidth',0.1,'DisplayName',sprintf('Decision variable for other states, precision: %.1f',precis));
        end
        plot(squeeze(mean(true_states_evidence,3)),'Color',precis_colors(color_iter,:),'LineWidth',2,'DisplayName',sprintf('Average trajectory for true state, precision: %.1f',precis));
        plot(squeeze(mean(sum(other_states_evidence,1),3)),'Color',precis_colors(color_iter,:),'LineStyle','--','LineWidth',2,'DisplayName',sprintf('Average trajectory for other states, precision: %.1f',precis));
        
    else
        
        all_means_true(i,:) = squeeze(mean(true_states_evidence,3));
        all_sems_true(i,:) = squeeze(std(true_states_evidence,0,3))./sqrt(num_repeats);
        
        all_means_other(i,:) = squeeze(mean(sum(other_states_evidence,1),3));
        all_sems_other(i,:) = squeeze(std(sum(other_states_evidence,1),0,3))./sqrt(num_repeats);
        
    end
        
   
    clear MDP
    clear x
    color_iter = color_iter + 1;
end

if num_repeats > 10
    
    lineProps.col = mat2cell(precis_colors,ones(size(precis_colors,1),1),3);
    
    mseb(1:T,all_means_true,all_sems_true,lineProps)
    
    lineProps.style = '--';
    hold on;
    mseb(1:T,all_means_other,all_sems_other,lineProps);
    
    legend([precis_labels_true,precis_labels_other])
    
end

xlim([1,T]);
legend('show')
xlabel('Sample index')
ylabel('Decision variable (x)')
