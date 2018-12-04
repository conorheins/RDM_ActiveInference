% step through single trials and print out posterior beliefs about outcomes and hidden
% states

num_epochs = size(MDPresult,2);

outcome_list = {'null','UP','RIGHT','DOWN','LEFT'};

for epoch = 1:num_epochs
    
    fprintf('Trial: %d ================================ \n',epoch)
    MDP = MDPresult(epoch);
    
    view_trial_trajectory_new(MDPresult,epoch,all_configs)
    
    
    for t = 1:MDP.T
        
        fprintf('Trial: %d, Timestep: %d\n',epoch,t);
        
        if MDP.o(3,t) >= 6
            switch MDP.o(3,t)
                case 6
                    fprintf('CHOICE: UP_RIGHT\n')
                case 7
                    fprintf('CHOICE: RIGHT_DOWN\n')
                case 8
                    fprintf('CHOICE: DOWN_LEFT\n')
                case 9
                    fprintf('CHOICE: LEFT_UP\n')
            end
        else
            
            outcome_beliefs = MDP.O{1,t};
            [sorted_beliefs,sort_idx] = sort(outcome_beliefs,'descend');
            
            fprintf('Sampled outcome: %s\n',outcome_list{MDP.o(1,t)})
            
            if sort_idx(1) > 1
                fprintf('Believed outcomes:\n')
                fprintf('%s: %.2f percent, %s: %.2f percent, %s: %.2f percent, %s: %.2f percent\n',outcome_list{sort_idx(1)},sorted_beliefs(1)*100,...
                    outcome_list{sort_idx(2)},sorted_beliefs(2)*100,outcome_list{sort_idx(3)},sorted_beliefs(3)*100,outcome_list{sort_idx(4)},sorted_beliefs(4)*100);
            else
                fprintf('Null\n')
            end
            
        end
                
        UR_beliefs = sum(MDP.xn{1}(end,scene_idx{1},t,t),2);
        RD_beliefs = sum(MDP.xn{1}(end,scene_idx{2},t,t),2);
        DL_beliefs = sum(MDP.xn{1}(end,scene_idx{3},t,t),2);
        LU_beliefs = sum(MDP.xn{1}(end,scene_idx{4},t,t),2);
        
                
        fprintf('Believed states:\n')
        fprintf('UP_RIGHT: %.2f, RIGHT_DOWN: %.2f, DOWN_LEFT, %.2f, LEFT_UP, %.2f\n',UR_beliefs*100,RD_beliefs*100,DL_beliefs*100,LU_beliefs*100);
        fprintf('=========================================\n')
                
        pause;
        
    end
    
    close gcf;
    
end
            
            
            
        
        
