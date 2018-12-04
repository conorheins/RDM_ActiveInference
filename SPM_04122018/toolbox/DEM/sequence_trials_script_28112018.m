chi_values = exp(linspace(log(0.01),log(1),25)); % negative threshold on the entropy of posterior beliefs (uncertainty in hierarchical schemes)
noise_values = exp(linspace(log(0.05),log(4),11)); % noise values for parameterizing likelihood precision
prior_precision = 0.25; % prior beliefs
N = 25; % number of trials per chi-value/precision-value pair

[ all_configs,scene_idx ] = generate_scenes();

for chi_i = 1:length(chi_values)

    for p_i = 1:length(noise_values)
            
        % randomly initialize agents' beliefs
        prior_scene_beliefs = cell(N,1);
        idx_tmp = randi([1 4],N,1);
        prior_scene_beliefs(idx_tmp == 1) = {'UP_RIGHT'};
        prior_scene_beliefs(idx_tmp == 2) = {'RIGHT_DOWN'};
        prior_scene_beliefs(idx_tmp == 3) = {'DOWN_LEFT'};
        prior_scene_beliefs(idx_tmp == 4) = {'LEFT_UP'};
        clear idx_temp;
        
        % randomly initialize true hidden states
        true_scenes = randi([1 48],N,1);
        
        for trial = 1:N
            
            MDPdeep = initialize_MDPdeep(prior_scene_beliefs{trial},prior_precision,all_configs,scene_idx,true_scenes(trial));
            MDPshallow = initialize_MDPshallow_witha_newp(noise_values(p_i),32);
            MDPshallow.chi = chi_values(chi_i);
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
        
        evidence_accum_curves = zeros(N,4,32);
        posterior_belief_curves = zeros(N,4,32);
        evidence_true = zeros(N,32);
        posterior_true = zeros(N,32);
        true_patterns = zeros(N,1);
        
        for trial = 1:N
            
            first_quadrant_visit = find(MDPresult(trial).o(1,:)>1,1);

            mdp_visit = MDPresult(trial).mdp(first_quadrant_visit);
            true_patterns(trial) = mdp_visit.s(1)-1;
            
            % evidence accumulation using log-Bayes-factor style
            % accumulation
            A = mdp_visit.A{1}(2:end,2:end);
            evidence_vecs = -ones(4) + 2*eye(4); 
            x = zeros(4,32);
            for s_i = 1:4
                x(s_i,1) = evidence_vecs(s_i,:)*log(A(mdp_visit.o(1)-1,:)');
            end
            
            if length(mdp_visit.o)>1
                for s_i = 1:4
                    for t = 2:length(mdp_visit.o)
                        x(s_i,t) = x(s_i,t-1) + evidence_vecs(s_i,:) * log(A(mdp_visit.o(t)-1,:)');
                    end
                end
            end
            
            evidence_accum_curves(trial,:,:) = x;
            evidence_true(trial,:) = x(true_patterns(trial),:);
            
            % evidence accumulation via accumulation of posterior beliefs
            % over time
            x = zeros(4,32);

            for s_i = 1:4
                xn_temp = diag(squeeze(mdp_visit.xn{1}(end,s_i+1,:,:)));
                x(s_i,1:length(xn_temp)) = cumsum(xn_temp);
            end
            
            posterior_belief_curves(trial,:,:) = x;
            posterior_true(trial,:) = x(true_patterns(trial),:);
            
        end
     
         save(sprintf('C:/RDPsearch_Hierarch9_results/RDPsearch_Hierarch9_ChiId%d_PId%d.mat',chi_i,p_i),'MDPresult','chi_values','noise_values',...
             'all_configs','scene_idx','prior_scene_beliefs','true_scenes','evidence_accum_curves','evidence_true','posterior_belief_curves','posterior_true','-v7.3')

         clear MDPresult;

    end
end
