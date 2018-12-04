function dotParams = create_dotParams(numPatterns,coh,speed,direction,dotSize)
% CREATE_DOTPARAMS Creates parameter struct for dot patterns, used in
%                   subsequent 'displayDots' function
% INPUTS: numPatterns -- number of different random dot patterns to create
%                        parameters for
%         coh         -- coherence of motion, either a single number
%                        (assumed to extend to all patterns) or a vector
%                        with length == numPatterns
%         speed       -- speed of motion, either a single number of vector,
%                        as with coh above^
%         direction   -- direction of motion, single number or vector as
%                        above
%         dotSize     -- size of dots in pattern, either single number or
%                        vector as above

if numPatterns == 1
    
    dotParams.coh = coh(1);
    dotParams.speed = speed(1);
    dotParams.direction = direction(1);
    dotParams.dotSize = dotSize(1);
    
elseif numPatterns > 1
    
    for patt_i = 1:numPatterns
        if length(coh) == 1
            dotParams(patt_i).coh = coh;
        elseif length(coh) == numPatterns
            dotParams(patt_i).coh = coh(patt_i);
        else
            warning('Coherence parameters not formatted correctly, just using first value in coherence vector')
            dotParams(patt_i).coh = coh(1);
        end
        
        if length(speed) == 1
            dotParams(patt_i).speed = speed;
        elseif length(speed) == numPatterns
            dotParams(patt_i).speed = speed(patt_i);
        else
            warning('Speed parameters not formatted correctly, just using first value in speed vector')
            dotParams(patt_i).speed = speed(1);
        end
        
        if length(direction) == 1
            dotParams(patt_i).direction = direction;
        elseif length(direction) == numPatterns
            dotParams(patt_i).direction = direction(patt_i);
        else
            warning('Direction parameters not formatted correctly, just using first value in direction vector')
            dotParams(patt_i).direction = direction(1);
        end
        
        if length(dotSize) == 1
            dotParams(patt_i).dotSize = dotSize;
        elseif length(dotSize) == numPatterns
            dotParams(patt_i).direction = direction(patt_i);
        else
            warning('dotSize parameters not formatted correctly, just using first value in dotSize vector');
        end
        
    end
    
end




    
end
        