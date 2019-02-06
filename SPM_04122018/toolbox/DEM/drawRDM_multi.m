function [ movie ] = drawRDM_multi(dotParams,T)
%DRAWRDM_MULTI Function that animates multiple random dot motion (RDM)
%patterns using parameters stored in structure array dotParams, which has
%one row per RDM
%   INPUTS: 1. dotParams, structure array containing the following fields (per
%           row/RDM pattern)
%               -center: coordinates of the center [x y] of the corresponding RDM
%               pattern
%               -apSize: size in [w h] of the aperture
%               -edge_spillover: how much the quadrant should respectively extend beyond the
%               width and height of the aperture, in pixels
%               -nDots: number of dots in the RDM
%               -speed: speed of motion in pixels/second
%               -direction: direction (0 to 360 degrees) of motion
%               -coherence: coherence of motion (% of dots moving in direction)
%           2. T: duration of movie, in frames



%%
nPatterns = length(dotParams);

quadrantsX = zeros(4,nPatterns);
quadrantsY = zeros(4,nPatterns);

for patt_i = 1:nPatterns
    
    quadrantsX(:,patt_i) = [dotParams(patt_i).center(1) - dotParams(patt_i).apSize(1) - dotParams(patt_i).edge_spillover(1);
         dotParams(patt_i).center(1) - dotParams(patt_i).apSize(1) - dotParams(patt_i).edge_spillover(1);
         dotParams(patt_i).center(1) + dotParams(patt_i).apSize(1) + dotParams(patt_i).edge_spillover(1);
         dotParams(patt_i).center(1) + dotParams(patt_i).apSize(1) + dotParams(patt_i).edge_spillover(1)];
     
    quadrantsY(:,patt_i) = [dotParams(patt_i).center(2) - dotParams(patt_i).apSize(2) - dotParams(patt_i).edge_spillover(2);
         dotParams(patt_i).center(2) + dotParams(patt_i).apSize(2) + dotParams(patt_i).edge_spillover(2);
         dotParams(patt_i).center(2) + dotParams(patt_i).apSize(2) + dotParams(patt_i).edge_spillover(2);
         dotParams(patt_i).center(2) - dotParams(patt_i).apSize(2) - dotParams(patt_i).edge_spillover(2)];

end


%%
x = [];
y = [];

dx = [];
dy = [];

l = [];
r = [];
b = [];
t = [];

apSizes = [];

lifetimes = [];

centers = [];

for patt_i = 1:nPatterns
    
    x = [x ,...
        (rand(1,dotParams(patt_i).nDots)-.5) * dotParams(patt_i).apSize(1) + dotParams(patt_i).center(1)];
    
    y = [y ,...
        (rand(1,dotParams(patt_i).nDots)-.5) * dotParams(patt_i).apSize(2) + dotParams(patt_i).center(2)];
    
    
    % get indices of dots that will be moving coherently
    coh_idx = rand(1,dotParams(patt_i).nDots) < dotParams(patt_i).coherence/100;
    
    % initialize motion vectors of all dots to random directions
    
%     rand_directions = 90.*randi(4,1,dotParams(patt_i).nDots) .* pi/180;
    rand_directions = 360.*rand(1,dotParams(patt_i).nDots) .* pi/180;
    
    dx_temp = 0.2*sin(rand_directions);
    dy_temp = -0.2*cos(rand_directions);
    
    dx_temp(coh_idx) = 0.2*sin(dotParams(patt_i).direction * pi/180);
    dy_temp(coh_idx) = -0.2*cos(dotParams(patt_i).direction * pi/180);
    
    dx = [dx,dx_temp];
    dy = [dy,dy_temp];
    
%     dx = [dx, 0.2*sin(dotParams(patt_i).direction * pi/180) * ones(1,dotParams(patt_i).nDots)];
%     dy = [dy, -0.2*cos(dotParams(patt_i).direction * pi/180) * ones(1,dotParams(patt_i).nDots)]; 
    
    l = [l, dotParams(patt_i).center(1)-dotParams(patt_i).apSize(1)/2 * ones(1,dotParams(patt_i).nDots)];
    r = [r, dotParams(patt_i).center(1)+dotParams(patt_i).apSize(1)/2 * ones(1,dotParams(patt_i).nDots)];
    b = [b, dotParams(patt_i).center(2)-dotParams(patt_i).apSize(2)/2 * ones(1,dotParams(patt_i).nDots)];
    t = [t, dotParams(patt_i).center(2)+dotParams(patt_i).apSize(2)/2 * ones(1,dotParams(patt_i).nDots)];
    
    apSizes = [apSizes;
               repmat(dotParams(patt_i).apSize,dotParams(patt_i).nDots,1)];
           
    lifetimes = [lifetimes, dotParams(patt_i).lifetime*ones(1,dotParams(patt_i).nDots)];
    
    centers = [centers;
                repmat(dotParams(patt_i).center,dotParams(patt_i).nDots,1)];
   
end
  
%%
       
lives =    ceil(rand(1,length(lifetimes)).*lifetimes);

movie = [];

for timestep = 1:T
    
    patch(quadrantsX,quadrantsY,'black')    

    x = x + dx;
    y = y + dy;
        
    x(x<l) = x(x<l) + apSizes(x<l,1)';
    x(x>r) = x(x>r) - apSizes(x>r,1)';
    y(y<b) = y(y<b) + apSizes(y<b,2)';
    y(y>t) = y(y>t) - apSizes(y>t,2)';
    
    % increment the lives of the dots
    lives = lives + 1;
    
    %find the 'dead' dots
    deadDots = mod(lives,lifetimes)==0;
     
    %replace the positions of the dead dots to a random location
    x(deadDots) = (rand(1,sum(deadDots))-.5).*apSizes(deadDots,1)' + centers(deadDots,1)';
    y(deadDots) = (rand(1,sum(deadDots))-.5).*apSizes(deadDots,2)' + centers(deadDots,2)';  
    
    hold on; scatter(x,y,'ko','MarkerFaceColor','w');

    curr_frame = getframe(gcf);
    
    movie = cat(4,movie,curr_frame.cdata);
    
    close gcf;    

end

