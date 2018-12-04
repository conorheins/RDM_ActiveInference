

[fnam,fdir] = uigetfile('*.mat');

load(fullfile(fdir,fnam));

trial_idx = 1;

M = MDP_beliefs_video(MDPresult,trial_idx,all_configs);

% video_folder = '/Users/conorheins/Documents/MATLAB/spm12_r6906_mexmaci64/toolbox/DEM/';
% video_name = sprintf('RDPsearch_trial%d_%s.gif',trial_idx,datestr(now,'ddmmyy'));

video_folder = '/Users/conorheins/Documents/Presentations/';
video_name = 'ExampleGifIncorrect5_08Uncertainty_06Prior.gif';

for i = 1:length(M)
    [img,cmap] = rgb2ind(M(i).cdata,256); 
    if i == 1 
        imwrite(img,cmap,fullfile(video_folder,video_name),'gif','LoopCount',Inf,'DelayTime',0.2); 
    elseif i == length(M)
        for repeat = 1:10 % hold the last frame for a bit
            imwrite(img,cmap,fullfile(video_folder,video_name),'gif','WriteMode','append','DelayTime',0.5);
        end
    else
        imwrite(img,cmap,fullfile(video_folder,video_name),'gif','WriteMode','append','DelayTime',0.2);
    end
end
