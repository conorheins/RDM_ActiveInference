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
all_errors_true = zeros(length(precisions2test),T);

all_means_other = zeros(length(precisions2test),T);
all_errors_other = zeros(length(precisions2test),T);

for i = 1:length(precisions2test)
    precis = precisions2test(i);
    D{1} = ones(num_outcomes,1);
    A{1} = spm_softmax(exp(precis)*log(eye(num_outcomes)+4));
    B{1}(:,:,1) = eye(num_outcomes);           % Actions do nothing (the two B matrices here are only included to prevent the scheme entering 'HMM mode')
    B{1}(:,:,2) = eye(num_outcomes);
    C{1} = zeros(num_outcomes,1);
    
    for rep_i = 1:num_repeats
        MDP(1,rep_i).A = A;
        MDP(1,rep_i).B = B;
        MDP(1,rep_i).C = C;
        MDP(1,rep_i).D = D;
        MDP(1,rep_i).T = T;
        MDP(1,rep_i).s = randi(num_outcomes)*ones(1,T);
    end

    mdp = spm_MDP_check(MDP);
    MDP = spm_MDP_VB_X(mdp);
    
    x = zeros(num_outcomes,T,num_repeats);
%     all_beliefs = zeros(size(x));

    evidence_vecs = -ones(num_outcomes)/(num_outcomes-1) + ((num_outcomes)/(num_outcomes-1))*eye(num_outcomes);

    for rep_i = 1:num_repeats
        
        for s_i = 1:num_outcomes
            x(s_i,1,rep_i) = evidence_vecs(s_i,:)*log(MDP(1,rep_i).A{1}(MDP(1,rep_i).o(1),:)');
        end
        
        for s_i = 1:num_outcomes
            for t = 2:T
                x(s_i,t,rep_i) = x(s_i,t-1) + evidence_vecs(s_i,:) * log(MDP(1,rep_i).A{1}(MDP(1,rep_i).o(t),:)');
            end
        end
    
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
    
    true_states = zeros(1,T,num_repeats);
    other_states = zeros(1,T,num_repeats);
    for rep_i = 1:num_repeats
        for s_i = 1:num_outcomes
            
            if s_i == MDP(1,rep_i).s(1)
                true_states(1,:,rep_i) = x(s_i,:,rep_i);
            else
                other_states(1,:,rep_i) = x(s_i,:,rep_i);
            end
            
        end
    end
    
   
    if num_repeats <= 10
        for rep_i = 1:num_repeats
            plot(squeeze(true_states(1,:,rep_i)),'Color',precis_colors(color_iter,:),'LineWidth',0.1,'DisplayName',sprintf('Decision variable for true state, precision: %.1f',precis));
            hold on;
            plot(squeeze(other_states(1,:,rep_i)),'Color',precis_colors(color_iter,:),'LineWidth',0.1,'DisplayName',sprintf('Decision variable for other states, precision: %.1f',precis));
        end
        plot(squeeze(mean(true_states,3)),'Color',precis_colors(color_iter,:),'LineWidth',2,'DisplayName',sprintf('Average trajectory for true state, precision: %.1f',precis));
        plot(squeeze(mean(other_states,3)),'Color',precis_colors(color_iter,:),'LineStyle','--','LineWidth',2,'DisplayName',sprintf('Average trajectory for other states, precision: %.1f',precis));
        
    else
        
        all_means_true(i,:) = squeeze(mean(true_states,3));
        all_errors_true(i,:) = squeeze(std(true_states,0,3))./sqrt(num_repeats);
        
        all_means_other(i,:) = squeeze(mean(other_states,3));
        all_errors_other(i,:) = squeeze(std(other_states,0,3))./sqrt(num_repeats);
        
    end
        
   
    clear MDP
    clear x
    color_iter = color_iter + 1;
end

if num_repeats > 10
    
    lineProps.col = mat2cell(precis_colors,ones(size(precis_colors,1),1),3);
    
    mseb(1:T,all_means_true,all_errors_true,lineProps)
    
    lineProps.style = '--';
    hold on;
    mseb(1:T,all_means_other,all_errors_other,lineProps);
    
    legend([precis_labels_true,precis_labels_other])
    
end

xlim([1,T]);
legend('show')
xlabel('Sample index')
ylabel('Decision variable (x)')
