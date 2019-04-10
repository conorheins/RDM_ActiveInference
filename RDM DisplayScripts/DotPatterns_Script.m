% Clear the workspace and the screen
sca;
close all;
clearvars;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer.
screens = Screen('Screens'); 

% sub_rect = [0 0 500 500];
sub_rect = []; % use full screen

monWidth = 50; % horizontal width of monitor being projected to
viewDist = 50; % distance of subject from monitor
% [ ScreenInfo ] = ScreenSetUp(screens(2),monWidth,viewDist,sub_rect); % use the connected external monitor
[ ScreenInfo ] = ScreenSetUp(screens(1),monWidth,viewDist,sub_rect); % use the connected external monitor


numPatterns = 2;
coh = [0.5 0.5];
speed = 1;
direction = [90 0];
dotSize = 3.5;

% numPatterns = 1;
% coh = 1;
% speed = 1;
% direction = 90;
% dotSize = 3.5;

dotParams = create_dotParams(numPatterns,coh,speed,direction,dotSize);

M = displayDots(ScreenInfo,dotParams,5,0);
