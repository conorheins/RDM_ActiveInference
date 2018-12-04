AssertOpenGL;
duration = 5; % how long to show the dots in seconds

curScreen = 0;
dontclear = 0;

[curWindow, screenRect] = Screen('OpenWindow', curScreen, [0,0,0],[],32, 2);

% Enable alpha blending with proper blend-function. We need it for drawing
% of smoothed points.
Screen('BlendFunction', curWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

fix_center = [screenRect(3), screenRect(4)]/2;

centers = [ [screenRect(3)/4, screenRect(4)/2] ; [(3*screenRect(3))/4, screenRect(4)/2] ];

Screen('Flip',curWindow,0,dontclear);

spf =Screen('GetFlipInterval', curWindow); % seconds per frame
monRefresh = 1/spf; % frames per second

% mon_horizontal_cm  	= 14.9865; % cm for window of these dimensions [300 50 1000 700]
mon_horizontal_cm   = 35;
view_dist_cm 		= 50;

% MAKE SURE APD IS USED CORRECTLY EVERYWHERE!!!!
% diameter/length of side of aperture
apD = 50;

% Everything is initially in coordinates of visual degrees, convert to pixels
% (pix/screen) * (screen/rad) * rad/deg
ppd = pi * screenRect(3) / atan(mon_horizontal_cm/view_dist_cm/2)  / 360;

d_ppd = floor(apD/10 * ppd);

% Dot stuff
coh = 0.7;
speed = 5;
direction = 90;
dotSize = 5;

maxDotsPerFrame = 150; % By trial and error and depends on graphics card

% ndots is the number of dots shown per video frame. Dots are placed in a
% square of the size of the aperture.
%   Size of aperture = Apd*Apd/100  sq deg
%   Number of dots per video frame = 16.7 dots per sq.deg/sec,
% When rounding up, do not exceed the number of dots that can be plotted in
% a video frame.
ndots = min(maxDotsPerFrame, ceil(16.7 * apD .* apD * 0.01 / monRefresh));

% dxdy is an N x 2 matrix that gives jumpsize in units on 0..1
%   deg/sec * Ap-unit/deg * sec/jump = unit/jump
dxdy = repmat(speed * 10/apD * (3/monRefresh) ...
    * [cos(pi*direction/180.0) -sin(pi*direction/180.0)], ndots,1);

% ARRAYS, INDICES for loop
ss1 = rand(ndots*3, 2); % array of dot positions raw [xposition, yposition] for first dot pattern
ss2 = rand(ndots*3, 2); % array of dot positions raw [xposition, yposition] for second dot pattern

% Divide dots into three sets
Ls = cumsum(ones(ndots,3)) + repmat([0 ndots ndots*2], ndots, 1);
loopi = 1; % Loops through the three sets of dots

% Show for how many frames
continue_show = round(duration*monRefresh);

priorityLevel = MaxPriority(curWindow,'KbCheck');
Priority(priorityLevel);

% THE MAIN LOOP

frames = 0;

while continue_show
    % Get ss & xs from the big matrices. xs and ss are matrices that have
    % stuff for dots from the last 2 positions + current.
    % Ls picks out the previous set (1:5, 6:10, or 11:15)
    Lthis  = Ls(:,loopi); % Lthis picks out the loop from 3 times ago, which
    % is what is then moved in the current loop
    this_s1 = ss1(Lthis,:);  % this is a matrix of random #s - starting positions for the first dot pattern
    this_s2 = ss2(Lthis,:);  % this is a matrix of random #s - starting positions for the second dot pattern
    
    % THE FOLLOWING IS DONE SEPARATELY FOR EACH DOT PATTERN
    % 1 group of dots are shown in the first frame, a second group are shown
    % in the second frame, a third group shown in the third frame. Then in
    % the next frame, some percentage of the dots from the first frame are
    % replotted according to the speed/direction and coherence, the next
    % frame the same is done for the second group, etc.
    
    % Update the loop pointer
    loopi = loopi+1;
    
    if loopi == 4
        loopi = 1;
    end
    
    % Compute new locations
    % L are the dots that will be moved
    L1 = rand(ndots,1) < coh;
    L2 = rand(ndots,1) < coh;
    this_s1(L1,:) = this_s1(L1,:) + dxdy(L1,:);	% Offset the selected dots for the first dot pattern
    this_s2(L2,:) = this_s1(L2,:) + dxdy(L2,:);    % Offset the selected dots for the second dot pattern
    
    if sum(~L1) > 0  % if not 100% coherence
        this_s1(~L1,:) = rand(sum(~L1),2);	% get new random locations for the rest
    end
    
    if sum(~L2) > 0 
        this_s1(~L2,:) = rand(sum(~L2),2); % get new random locations for the rest
    end
    
    % Wrap around, Dot Pattern 1 - check to see if any positions are greater than one or
    % less than zero which is out of the aperture, and then replace with a
    % dot along one of the edges opposite from direction of motion.
    
    N1 = sum((this_s1 > 1 | this_s1 < 0)')' ~= 0;
    if sum(N1) > 0
        xdir = sin(pi*direction/180.0);
        ydir = cos(pi*direction/180.0);
        % Flip a weighted coin to see which edge to put the replaced dots
        if rand < abs(xdir)/(abs(xdir) + abs(ydir))
            this_s1(find(N1==1),:) = [rand(sum(N1),1) (xdir > 0)*ones(sum(N1),1)];
        else
            this_s1(find(N1==1),:) = [(ydir < 0)*ones(sum(N1),1) rand(sum(N1),1)];
        end
    end
    
    % Wrap around, Dot Pattern 2 - check to see if any positions are greater than one or
    % less than zero which is out of the aperture, and then replace with a
    % dot along one of the edges opposite from direction of motion.
    
    N2 = sum((this_s2 > 1 | this_s2 < 0)')' ~= 0;
    if sum(N2) > 0
        xdir = sin(pi*direction/180.0);
        ydir = cos(pi*direction/180.0);
        % Flip a weighted coin to see which edge to put the replaced dots
        if rand < abs(xdir)/(abs(xdir) + abs(ydir))
            this_s2(find(N2==1),:) = [rand(sum(N2),1) (xdir > 0)*ones(sum(N2),1)];
        else
            this_s2(find(N2==1),:) = [(ydir < 0)*ones(sum(N2),1) rand(sum(N2),1)];
        end
    end
    
    % Convert to stuff we can actually plot
    this_x1(:,1:2) = floor(d_ppd(1) * this_s1); % pix/ApUnit for first dot pattern
    this_x2(:,1:2) = floor(d_ppd(1) * this_s2); % pix/ApUnit for second dot pattern
    
    
    % This assumes that zero is at the top left, but we want it to be in the
    % center, so shift the dots up and left, which just means adding half of
    % the aperture size to both the x and y direction.
    dot_show1 = (this_x1(:,1:2) - d_ppd/2)';
    dot_show2 = (this_x2(:,1:2) - d_ppd/2)';
    
    % After all computations, flip
    Screen('Flip', curWindow,0,dontclear);
    % Now do next drawing commands
    
    Screen('DrawDots', curWindow, dot_show1, dotSize, [255 255 255], centers(1,:));
    Screen('DrawDots', curWindow, dot_show2,dotSize, [255 255 255], centers(2,:));
    
    Screen('DrawDots', curWindow, [0; 0], 10, [255 0 0], fix_center, 1);
    
    % Presentation
    Screen('DrawingFinished',curWindow,dontclear);
    
    frames = frames + 1;
    
    if frames == 1
        start_time = GetSecs;
    end
    
    % Update the arrays so xor works next time
    xs1(Lthis, :) = this_x1;
    xs2(Lthis, :) = this_x2;
    ss1(Lthis, :) = this_s1;
    ss2(Lthis, :) = this_s2;
    
    % Check for end of loop
    continue_show = continue_show - 1;
    
end

% Present last dots
Screen('Flip', curWindow,0,dontclear);

% Erase last dots
Screen('DrawingFinished',curWindow,dontclear);
Screen('Flip', curWindow,0,dontclear);

Screen('CloseAll'); % Close display windows
Priority(0); % Shutdown realtime mode.
ShowCursor; % Show cursor again, if it has been disabled.