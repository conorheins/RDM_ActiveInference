function [ ScreenInfo ] = ScreenSetUp(curScreen,monWidth,viewDist,sub_rect)
% ScreenSetUp Function that initializes experimental screen, collects info
% for future use.
%   INPUTS: curScreen -- index of the screen to create display upon
%           monWidth  -- width of the screen, in cm. Will be changed if
%                        sub_rect is not the whole screen
%           viewDist --  viewing distance from experimental subject to
%                        screen, in cm
%           sub_rect -- [1 x 4] vector determining size of sub-window to
%                       draw within curScreen 

AssertOpenGL;

dontclear = 0;

[curWindow, screenRect] = Screen('OpenWindow', curScreen, [0,0,0],[],32, 2);
if nargin < 4 || isempty(sub_rect)
    ScreenInfo.screenRect = screenRect;
elseif exist('sub_rect','var') && sum(sub_rect == screenRect) < 4
    sca;
    Screen('Preference', 'SkipSyncTests', 1);
    [curWindow, screenRect_sub] = Screen('OpenWindow',curScreen,[0,0,0],sub_rect, 32, 2);
    monWidth = monWidth * (screenRect_sub(3) - screenRect_sub(1))/(screenRect(3) - screenRect(1));
    ScreenInfo.screenRect = screenRect_sub;
end

% Enable alpha blending with proper blend-function. We need it for drawing
% of smoothed points.
Screen('BlendFunction', curWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

Screen('Flip',curWindow,0,dontclear);

spf =Screen('GetFlipInterval', curWindow); % seconds per frame
ScreenInfo.monRefresh = 1/spf; % frames per second; 

% MAKE SURE APD IS USED CORRECTLY EVERYWHERE!!!!
% diameter/length of side of aperture
ScreenInfo.apD = 50;

% Everything is initially in coordinates of visual degrees, convert to pixels
% (pix/screen) * (screen/rad) * rad/deg
ScreenInfo.ppd = pi * ScreenInfo.screenRect(3) / atan(monWidth/viewDist/2)  / 360;

ScreenInfo.d_ppd = floor(ScreenInfo.apD/10 * ScreenInfo.ppd);

ScreenInfo.curWindow = curWindow;

end

