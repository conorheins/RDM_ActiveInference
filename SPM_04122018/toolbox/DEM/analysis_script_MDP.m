%% analyze MDP behavior & beliefs script

[fnam,fdir] = uigetfile('*.mat'); % choose data file to analyze

load(fullfile(fdir,fnam));

num_trials = length(MDP);

[~,Scenes] = initialize_searchMDP_4(0.1,1,'UP_RIGHT',0.25); % just use this line to get the scene identities

for t = 1:num_trials
    view_trial_trajectory_new(MDP,t,Scenes)    
    title(sprintf('Trial %d, Incoherence Level: %.2f, Prior Belief: %.2f in %s\n',t,noise_val(t),prior_scene_prob(t),strrep(prior_scenes{t},'_',' ')));
    pause;
    close gcf;
end

behav_vars = {'Accuracy','Saccades','RT'};

stats_table = zeros(num_trials,length(behav_vars));

for t = 1:num_trials
    o = MDP(t).o(1,:);
    stats_table(t,strcmp(behav_vars,'Accuracy')) = double(any(o == 6) & ~any(o == 7));  % accuracy
    stats_table(t,strcmp(behav_vars,'Saccades')) = find([(o > 5), 1],1) - 1;            % number of saccades
    stats_table(t,strcmp(behav_vars,'RT')) = mean(MDP(t).rt);                           % reaction time
end

% relate frequency of revisits to quadrants to uncertainty and prior probability

True_Scenes = {'UP_RIGHT','RIGHT_DOWN','DOWN_LEFT','LEFT_UP'};
vars = {'Number Revisits','Uncertainty','True Scene','Believed Scene','Prior Probability of Belief','Agreement'};
stats_table = zeros(num_trials,length(vars));
for t = 1:num_trials
    
    % calculate repeated visits to quadrants
    o = MDP(t).o(1,:);
    for quadrant = 2:5
        if find(o > 5,1) == 2 & o(3) == o(2) % don't analyze trials where they just immediately saccaded to their prior belief and stuck there
            stats_table(t,strcmp(vars,'Number Revisits')) = NaN;
        elseif length(find(o == quadrant)) > 1
            stats_table(t,strcmp(vars,'Number Revisits')) = stats_table(t,strcmp(vars,'Number Revisits')) + length(find(o == quadrant));
        end
    end
    
    % fill out uncertainty values
    stats_table(t,strcmp(vars,'Uncertainty')) = noise_val(t);
    
    % fill out true scene IDs
    if sum(strcmp('U',Scenes(MDP(t).s(1,1),:))) == 1 && sum(strcmp('R',Scenes(MDP(t).s(1,1),:))) == 1
        stats_table(t,strcmp(vars,'True Scene')) = 1;
    elseif sum(strcmp('R',Scenes(MDP(t).s(1,1),:))) == 1 && sum(strcmp('D',Scenes(MDP(t).s(1,1),:))) == 1
        stats_table(t,strcmp(vars,'True Scene')) = 2;
    elseif sum(strcmp('D',Scenes(MDP(t).s(1,1),:))) == 1 && sum(strcmp('L',Scenes(MDP(t).s(1,1),:))) == 1
        stats_table(t,strcmp(vars,'True Scene')) = 3;
    elseif sum(strcmp('D',Scenes(MDP(t).s(1,1),:))) == 1 && sum(strcmp('L',Scenes(MDP(t).s(1,1),:))) == 1
        stats_table(t,strcmp(vars,'True Scene')) = 4;
    end
    
    % fill out believed scene IDs
    stats_table(t,strcmp(vars,'Believed Scene')) = find(strcmp(True_Scenes,prior_scenes(t)));
    
    % fill out prior belief probabilities
    stats_table(t,strcmp(vars,'Prior Probability of Belief')) = prior_scene_prob(t);
    
    stats_table(t,strcmp(vars,'Agreement')) = stats_table(t,strcmp(vars,'True Scene')) == stats_table(t,strcmp(vars,'Believed Scene'));
    
end

figure;
scatter3(stats_table(stats_table(:,strcmp(vars,'Agreement')) == 1,strcmp(vars,'Number Revisits')),stats_table(stats_table(:,strcmp(vars,'Agreement')) == 1,strcmp(vars,'Uncertainty')),...
    stats_table(stats_table(:,strcmp(vars,'Agreement')) == 1,strcmp(vars,'Prior Probability of Belief')),'bo')
hold on; scatter3(stats_table(stats_table(:,strcmp(vars,'Agreement')) == 0,strcmp(vars,'Number Revisits')),stats_table(stats_table(:,strcmp(vars,'Agreement')) == 0,strcmp(vars,'Uncertainty')),...
    stats_table(stats_table(:,strcmp(vars,'Agreement')) == 0,strcmp(vars,'Prior Probability of Belief')),'ro')        
xlabel('Number Revisits');
ylabel('Uncertainty')
zlabel('Prior Probability of Belief')
legend({'Correct prior','Incorrect prior'})
legend('show');


num_revisits_correct = stats_table(stats_table(:,strcmp(vars,'Agreement')) == 1 & ~isnan(stats_table(:,strcmp(vars,'Number Revisits'))), strcmp(vars,'Number Revisits'));
uncertainty_correct = noise_val( stats_table(:,strcmp(vars,'Agreement')) == 1 & ~isnan(stats_table(:,strcmp(vars,'Number Revisits'))));
prior_strength_correct = prior_scene_prob(stats_table(:,strcmp(vars,'Agreement')) == 1 & ~isnan(stats_table(:,strcmp(vars,'Number Revisits'))));

prior_thr = prctile(prior_strength_correct,50);
scatter(num_revisits_correct(prior_strength_correct < prior_thr),uncertainty_correct(prior_strength_correct < prior_thr),'b','DisplayName','Weak Priors')
hold on; scatter(num_revisits_correct(prior_strength_correct >= prior_thr),uncertainty_correct(prior_strength_correct >= prior_thr),'r','DisplayName','Strong Priors')
xlabel('Number Revisits')
ylabel('Uncertainty')
legend('show')


Pquants = quantile(prior_strength_correct,[0.333 0.666]);
Uquants = quantile(uncertainty_correct,[0.25 0.5 0.75]);

bar_graph_data = zeros(4,3,2);
for p_i = 1:(length(Pquants)+1)
    if p_i == 1
        cond1 = prior_strength_correct < Pquants(p_i);
    elseif p_i == (length(Pquants)+1)
        cond1 = prior_strength_correct > Pquants(p_i-1);
    else
        cond1 = prior_strength_correct > Pquants(p_i-1) & prior_strength_correct < Pquants(p_i); 
    end
    for u_i = 1:(length(Uquants)+1)
        if u_i == 1
            cond2 = uncertainty_correct < Uquants(u_i);
        elseif u_i == (length(Uquants)+1)
            cond2 = uncertainty_correct > Uquants(u_i-1);
        else
            cond2 = uncertainty_correct > Uquants(u_i - 1) & uncertainty_correct < Uquants(u_i);
        end
        bar_graph_data(u_i,p_i,1) = mean(num_revisits_correct(cond1 & cond2));
        bar_graph_data(u_i,p_i,2) = std(num_revisits_correct(cond1 & cond2))./sqrt(length(num_revisits_correct(cond1 & cond2)));
    end
end

% use barwitherr.m (File Exchange contribution) to plot error bars around
% grouped bar graph
barwitherr(bar_graph_data(:,:,2),bar_graph_data(:,:,1))

%% analyze systematic results

[fnam,fdir] = uigetfile('*.mat'); % choose data file to analyze

load(fullfile(fdir,fnam));

[~,Scenes] = initialize_searchMDP_4(0.1,1,'UP_RIGHT',0.25); % just use this line to get the scene identities
believed_scene = 'UP_RIGHT';

vars = {'Number Revisits','Accuracy','Saccades'};

stats_table = zeros([size(results),length(vars)]);

for ii = 1:length(noise_vals)
    for jj = 1:length(prior_prob)
        
        num_revisits = zeros(100,1);
        acc = zeros(100,1);
        num_saccades = zeros(100,1);
        for n = 1:100
            % calculate repeated visits to quadrants
            o = results{ii,jj,1}(n).o(1,:);
            for quadrant = 2:5
                if find(o > 5,1) == 2 & o(3) == o(2) % don't analyze trials where they just immediately saccaded to their prior belief and stuck there
                    num_revisits(n) = NaN;
                elseif length(find(o == quadrant)) > 1
                    num_revisits(n) = num_revisits(n) + length(find(o == quadrant));
                end
            end
            
            acc(n) = double(any(o == 6) & ~any(o == 7));           % accuracy
            num_saccades(n) = find([(o > 5), 1],1) - 1;            % number of saccades
            
        end
        
        stats_table(ii,jj,1,1) = mean(num_revisits(~isnan(num_revisits)));
        stats_table(ii,jj,1,2) = mean(acc);
        stats_table(ii,jj,1,3) = mean(num_saccades);
        
        num_revisits = zeros(100,1);
        acc = zeros(100,1);
        num_saccades = zeros(100,1);
        for n = 1:100
            % calculate repeated visits to quadrants
            o = results{ii,jj,2}(n).o(1,:);
            for quadrant = 2:5
                if find(o > 5,1) == 2 & o(3) == o(2) % don't analyze trials where they just immediately saccaded to their prior belief and stuck there
                    num_revisits(n) = NaN;
                elseif length(find(o == quadrant)) > 1
                    num_revisits(n) = num_revisits(n) + length(find(o == quadrant));
                end
            end
            
            acc(n) = double(any(o == 6) & ~any(o == 7));           % accuracy
            num_saccades(n) = find([(o > 5), 1],1) - 1;            % number of saccades
            
        end
        
        stats_table(ii,jj,2,1) = mean(num_revisits(~isnan(num_revisits)));
        stats_table(ii,jj,2,2) = mean(acc);
        stats_table(ii,jj,2,3) = mean(num_saccades);
        

    end
end

