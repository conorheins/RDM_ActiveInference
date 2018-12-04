function [ all_configs,scene_idx ] = generate_scenes()
% GENERATE_SCENES Generates all 48 scene configurations for the 'RDP (random dot pattern) Scene Construction' 
% study, as well as list of indices that map from 12 unique RDP-configurations to one scene identity. 
% The Scene names are as follows: 
% -'UP_RIGHT' (one RDP facing UP, one RDP facing RIGHT)
% -'RIGHT_DOWN' (one RDP facing RIGHT, one RDP facing DOWN)
% -'DOWN_LEFT' (one RDP facing DOWN, one RDP facing LEFT)
% -'LEFT_UP' (one RDP facing LEFT, one RDP facing UP)

num_scenes = 4; % four scenes
num_configs = 12; % 12 unique ways to populate four quadrants with the two RDPs comprising each scene

U_R = unique(perms([2 3 1 1]),'rows');
R_D = unique(perms([3 4 1 1]),'rows');
D_L = unique(perms([4 5 1 1]),'rows');
L_U = unique(perms([5 2 1 1]),'rows');

Scene = zeros(num_scenes*num_configs,size(U_R,2));
start_idx = [1:num_configs:(num_configs*num_scenes),num_configs*num_scenes];
scene_idx = cell(num_scenes,1);
scene_idx{1} = start_idx(1):(start_idx(2)-1);
scene_idx{2} = start_idx(2):(start_idx(3)-1);
scene_idx{3} = start_idx(3):(start_idx(4)-1);
scene_idx{4} = start_idx(4):start_idx(5);

for sc_i = 1:num_scenes
    Scene(scene_idx{1},:) = U_R;
    Scene(scene_idx{2},:) = R_D;
    Scene(scene_idx{3},:) = D_L;
    Scene(scene_idx{4},:) = L_U;
end

all_configs = num2cell(Scene);
all_configs([all_configs{:}] == 1) = {'n'}; 
all_configs([all_configs{:}] == 2) = {'U'};
all_configs([all_configs{:}] == 3) = {'R'};
all_configs([all_configs{:}] == 4) = {'D'};
all_configs([all_configs{:}] == 5) = {'L'};

end

