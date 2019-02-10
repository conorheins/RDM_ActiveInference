% draw RDM script, calls drawRDM_multi for animation

npatterns = 4;
T = 50;

centers = [-5 5 ; 
          -5 -5 ;
           5  5 ;
           5 -5];
       
apSizes = [4 4;
           4 4;
           4 4;
           4 4];
       
edge_spillovers = [0.5 0.5;
                    0.5 0.5;
                    0.5 0.5;
                    0.5 0.5;];
       

nDots = [10; 15; 20; 25];

speeds = [7; 7; 7; 7];

directions = [90; 180; 270; 360];

coherences = [100; 50; 20; 10];

lifetimes = [12; 12; 12; 12];

for patt_i = 1:npatterns
    
    dotParams(patt_i).center = centers(patt_i,:);
    dotParams(patt_i).apSize = apSizes(patt_i,:);
    dotParams(patt_i).edge_spillover = edge_spillovers(patt_i,:);
    dotParams(patt_i).nDots = nDots(patt_i);
    dotParams(patt_i).speed = speeds(patt_i);
    dotParams(patt_i).direction = directions(patt_i);
    dotParams(patt_i).coherence = coherences(patt_i);
    dotParams(patt_i).lifetime = lifetimes(patt_i);
    
end

movie = drawRDM_multi(dotParams,T);

for t = 1:T
    imshow(movie(:,:,:,t));
    pause(0.025);
end




