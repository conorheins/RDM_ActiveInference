function MDP_Drift_Diffusion_manyOutcomes(num_outcomes,precisions2test)

precis_colors = cool(length(precisions2test));

color_iter = 1;

for i = 1:length(precisions2test)
    precis = precisions2test(i);
    D{1} = ones(num_outcomes,1);
    A{1} = spm_softmax(exp(precis)*log(eye(num_outcomes)+4));
    B{1}(:,:,1) = eye(num_outcomes);           % Actions do nothing (the two B matrices here are only included to prevent the scheme entering 'HMM mode')
    B{1}(:,:,2) = eye(num_outcomes);
    C{1} = zeros(num_outcomes,1);
   
    MDP.A = A;
    MDP.B = B;
    MDP.C = C;
    MDP.D = D;
    MDP.T = 40;
    MDP.s = ones(1,MDP.T);
   
    mdp = spm_MDP_check(MDP);
    MDP = spm_MDP_VB_X(mdp);
    
    x = zeros(num_outcomes,MDP.T);
    evidence_vecs = -ones(num_outcomes) + 2*eye(num_outcomes);
    for s_i = 1:num_outcomes
        x(s_i,1) = evidence_vecs(s_i,:)*log(MDP.A{1}(MDP.o(1),:)');
    end
    
    for s_i = 1:num_outcomes
        for t = 2:MDP.T
            x(s_i,t) = x(s_i,t-1) + evidence_vecs(s_i,:) * log(MDP.A{1}(MDP.o(t),:)');
        end
    end
    
    all_beliefs = zeros(num_outcomes,MDP.T);
    for s_i = 1:num_outcomes
        all_beliefs(s_i,:) = cumsum(diag(squeeze(MDP.xn{1}(end,s_i,:,:))));
    end
    
    all_beliefs = spm_softmax(all_beliefs);
        
    negentropy = zeros(1,MDP.T);
    for t = 1:MDP.T
        negentropy(t) = all_beliefs(:,t)'*log(all_beliefs(:,t) + exp(-16));
    end
    
    plot(negentropy,'Color',precis_colors(color_iter,:),'DisplayName',sprintf('Negentropy of posterior for %d states, precision: %.1f',num_outcomes,precis));
    hold on;

%     for s_i = 1:numoutcomes
%         
%         if s_i == MDP.s(1)
%             plot(x(s_i,:),'Color',precis_colors(color_iter,:),'LineWidth',2,'DisplayName',sprintf('Evidence for true state, precision: %.1f',precis));
%         else
%             plot(x(s_i,:),'Color',precis_colors(color_iter,:),'LineWidth',0.5,'DisplayName',sprintf('Evidence for other states, precision: %.1f',precis));
%         end
%             
%         hold on;
%     end
    
    xlim([1 MDP.T])

    clear MDP
    clear x
    color_iter = color_iter + 1;
end

legend('show')
xlabel('Sample')
ylabel('Negentropy')
