function M = displayDots(ScreenInfo,dotParams,duration,make_movie_flag)
%displayDots Displays random dot pattern of given coherence, or a series of
%patches of random dot patterns, depending on whether dotParams is a
%structure or a struct array
%   INPUTS: -ScreenInfo: structure containing information about the display
%                        screen
%           -dotParams: structure or structure array containing
%                       patch-specific dot-pattern parameters (e.g. coherence,
%                       direction of motion, etc.) If dotParams is a struct
%                       array, it can be made of up to 4 patterns total
%                       (i.e. max length is 4)
%           -duration:  how long to show the dots in seconds
%           -make_movie_flag: flag for whether to save the output as a
%                       movie that can subsequently be written to disk as .gif, for
%                       instance

M = [];

numPatches = length(dotParams);
if numPatches > 4
    warning(sprintf('More than four dot patterns detected: only using first four entries of dotParams array\n'));
    dotParams = dotParams(1:4);
    numPatches = 4;
end

curWindow = ScreenInfo.curWindow;
screenRect = ScreenInfo.screenRect;
fix_center = [screenRect(3), screenRect(4)]/2;

centers = zeros(numPatches,2);
numRows = rem(numPatches,2) + floor(numPatches/2);
numCols = 2;
for i = 1:numPatches
    which_row = rem(i,2) + floor(i/2);
    which_col = abs(rem(i,2)-1) + 1;
    centers(i,:) = [which_col * (screenRect(3)/(numCols+1)),which_row*screenRect(4)/(numRows+1)];
end

dontclear = 0;

Screen('Flip',curWindow,0,dontclear);

% Dot stuff
coh = [dotParams(:).coh];
speed = [dotParams(:).speed];
direction = [dotParams(:).direction];
dotSize = [dotParams(:).dotSize];

apD = ScreenInfo.apD;
monRefresh = ScreenInfo.monRefresh;

maxDotsPerFrame = 20; % By trial and error and depends on graphics card

% ndots is the number of dots shown per video frame. Dots are placed in a
% square of the size of the aperture.
%   Size of aperture = Apd*Apd/100  sq deg
%   Number of dots per video frame = 16.7 dots per sq.deg/sec,
% When rounding up, do not exceed the number of dots that can be plotted in
% a video frame.

% ndots = min(maxDotsPerFrame, ceil(16.7 * apD .* apD * 0.01 / monRefresh));
ndots = 10;

% dxdy is an N x 2 matrix that gives jumpsize in units on 0..1
%   deg/sec * Ap-unit/deg * sec/jump = unit/jump

dxdy = zeros(ndots,2,numPatches);
ss = rand(ndots*3,2,numPatches);

for i = 1:numPatches
    
    dxdy(:,:,i) = repmat(speed(i) * 10/apD * (3/monRefresh) ...
        * [cos(pi*direction(i)/180.0), -sin(pi*direction(i)/180.0)],ndots,1);
    
end
    
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
    
    this_ss = ss(Lthis,:,:); % this is the array of random starting positions for the numPatches of RDPs
    

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
    
    all_dots_show = cell(1,numPatches);
    
    for i = 1:numPatches
        
        % Compute new locations
        % L are the dots that will be moved
        L = rand(ndots,1) < coh(i);
        this_ss(L,:,i) = this_ss(L,:,i) + dxdy(L,:,i); % Offset the selected dots for the ith dot pattern

        if sum(~L) > 0  % if not 100% coherence
            this_ss(~L,:,i) = rand(sum(~L),2);	% get new random locations for the rest
        end

        % Wrap around, Dot Pattern i - check to see if any positions are greater than one or
        % less than zero which is out of the aperture, and then replace with a
        % dot along one of the edges opposite from direction of motion.

        N = sum((this_ss(:,:,i) > 1 | this_ss(:,:,i) < 0)')' ~= 0;
        if sum(N) > 0
            xdir = sin(pi*direction(i)/180.0);
            ydir = cos(pi*direction(i)/180.0);
            % Flip a weighted coin to see which edge to put the replaced dots
            if rand < abs(xdir)/(abs(xdir) + abs(ydir))
                this_ss(N==1,:,i) = [rand(sum(N),1) (xdir > 0)*ones(sum(N),1)];
            else
                this_ss(N==1,:,i) = [(ydir < 0)*ones(sum(N),1) rand(sum(N),1)];
            end
        end
        % Convert to stuff we can actually plot
        this_x(:,1:2,i) = floor(ScreenInfo.d_ppd * this_ss(:,:,i)); % pix/ApUnit for first dot pattern

        dot_show = (this_x(:,1:2,i) - ScreenInfo.d_ppd/2)';
        
        all_dots_show{i} = dot_show;

%         % After all computations, flip
%         Screen('Flip', curWindow,0,dontclear);
%         % Now do next drawing commands
% 
%         Screen('DrawDots', curWindow, dot_show, dotSize(i), [255 255 255], centers(i,:));
%         Screen('DrawDots', curWindow, [0; 0], 10, [255 0 0], fix_center, 1); 
        
%         Screen('DrawingFinished',curWindow,dontclear);

    end
    
    % After all computations, flip
    Screen('Flip', curWindow,0,dontclear);
    % Now do next drawing commands
    
    for i = 1:numPatches
    
        Screen('DrawDots', curWindow, all_dots_show{i}, dotSize(i), [255 255 255], centers(i,:));
        Screen('DrawDots', curWindow, [0; 0], 10, [255 0 0], fix_center, 1);
        
    end
    
    % Presentation
    Screen('DrawingFinished',curWindow,dontclear);
    
    frames = frames + 1;
    
    if make_movie_flag
        imageArray = Screen('GetImage',curWindow,[300 300 1300 1300]);
        M = cat(3,M,imageArray);
    end
    
    if frames == 1
        start_time = GetSecs;
    end
    
    % Update the arrays so xor works next time
    xs(Lthis,:,:) = this_x;
    ss(Lthis,:,:) = this_ss;
    
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


end

