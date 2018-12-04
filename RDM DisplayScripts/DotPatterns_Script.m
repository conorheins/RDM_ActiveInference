% Clear the workspace and the screen
sca;
close all;
clearvars;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer.
screens = Screen('Screens'); 

monWidth = 50; % horizontal width of monitor being projected to
viewDist = 50; % distance of subject from monitor
[ ScreenInfo ] = ScreenSetUp(screens(2),monWidth,viewDist,[]); % use the connected external monitor


numPatterns = 2;
coh = [1 0.75];
speed = 2;
direction = [90 0];
dotSize = 7.5;

dotParams = create_dotParams(numPatterns,coh,speed,direction,dotSize);


M = displayDots(ScreenInfo,dotParams,5,1);
