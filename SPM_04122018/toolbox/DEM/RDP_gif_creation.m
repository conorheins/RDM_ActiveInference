

% [fnam,fdir] = uigetfile('*.mat');
% 
% load(fullfile(fdir,fnam));

prior_belief = 'RIGHT_DOWN';
prior_precis = 0.7;
true_scene = 3;
sens_precis = 4;

MDP = random_searchMDP(prior_belief,prior_precis,true_scene,sens_precis);

trial_idx = 1;

% M = MDP_beliefs_video(MDPresult,trial_idx,all_configs);
[M,dot_idx] = MDP_BeliefsVid_Hierarch(MDP,trial_idx);

% video_folder = '/Users/conorheins/Documents/MATLAB/spm12_r6906_mexmaci64/toolbox/DEM/';
% video_name = sprintf('RDPsearch_trial%d_%s.gif',trial_idx,datestr(now,'ddmmyy'));

video_folder = '/Users/conorheins/Documents/Presentations/';
% video_name = 'ExampleGifIncorrect5_08Uncertainty_06Prior.gif';
video_name = 'ExampleGifCorrect_400Precision_70Prior.gif';

for i = 1:length(M)
    [img,cmap] = rgb2ind(M(i).cdata,256); 
    if i == 1
        imwrite(img,cmap,fullfile(video_folder,video_name),'gif','LoopCount',Inf,'DelayTime',0.2); 
    elseif ismember(i,dot_idx)
        imwrite(img,cmap,fullfile(video_folder,video_name),'gif','WriteMode','append','DelayTime',0.1);
    elseif i == length(M)
        for repeat = 1:10 % hold the last frame for a bit
            imwrite(img,cmap,fullfile(video_folder,video_name),'gif','WriteMode','append','DelayTime',0.5);
        end
    else
        imwrite(img,cmap,fullfile(video_folder,video_name),'gif','WriteMode','append','DelayTime',0.2);
    end
end


% for i = 1:length(M)
%     [img,cmap] = rgb2ind(M(i).cdata,256); 
%     if i == 1 
%         imwrite(img,cmap,fullfile(video_folder,video_name),'gif','LoopCount',Inf,'DelayTime',0.2); 
%     elseif i == length(M)
%         for repeat = 1:10 % hold the last frame for a bit
%             imwrite(img,cmap,fullfile(video_folder,video_name),'gif','WriteMode','append','DelayTime',0.5);
%         end
%     else
%         imwrite(img,cmap,fullfile(video_folder,video_name),'gif','WriteMode','append','DelayTime',0.2);
%     end
% end
