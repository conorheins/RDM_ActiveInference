%% solve and evaluate performance for sequence of trials

% systematically test how uncertainty affects saccades under different
% degrees of prior confidence

noise_vals = [0.25 0.5 0.75 1.00 1.25 1.50];
prior_prob = [0.25 0.35 0.5 0.65 0.75 0.8];
believed_scene = 'UP_RIGHT';

results = cell(length(noise_vals),length(prior_prob),2);

for ii = 1:length(noise_vals)
    for jj = 1:length(prior_prob)
        
        for n = 1:100
            MDP_UR(n) = initialize_searchMDP_4(noise_vals(ii),randi([1 12],1,1),believed_scene,prior_prob(jj));
        end
        results{ii,jj,1} = spm_MDP_VB_X_distribO(MDP_UR);

        
        for n = 1:100
            MDP_RD(n) = initialize_searchMDP_4(noise_vals(ii),randi([13 24],1,1),believed_scene,prior_prob(jj));
        end
        results{ii,jj,2} = spm_MDP_VB_X_distribO(MDP_RD);

    end
end

save(sprintf('MDPresults_sys_%s.mat',datestr(now,'mmddyy')),'results','noise_vals','prior_prob','-v7.3')

