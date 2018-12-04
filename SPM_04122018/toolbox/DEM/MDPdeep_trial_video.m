function [ M ] = MDPdeep_trial_video(MDP,all_configs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% load images
%----------------------------------------------------------------------
load MDP_search_graphics_DOT
UP = double(imresize(UP,0.8));
RIGHT = double(imresize(RIGHT,0.8));
DOWN = double(imresize(DOWN,0.8));
LEFT = double(imresize(LEFT,0.8));
figure('Position',[400 300 600 450]);

% make a mask that represents currently-foveated outcome
null = zeros(size(UP)) + 1;
mask  = hamming(size(UP,1));
mask  = mask*mask';
for i = 1:3
    mask(:,:,i) = mask(:,:,1);
end

% locations
% x = x,y coordinates of 1) fixation point, 2) four quadrants, and 3) four
% categorization choices
%--------------------------------------------------------------------------
x = [0        0;
    -1       -1;
    -1        1;
    1        -1;
    1         1;
    -1.25   2.5;
    -0.42   2.5;
    0.42    2.5;
    1.25    2.5];

r = [-1,1]/2;

 
% plot cues
%--------------------------------------------------------------------------
Ni    = 1:size(MDP.xn{1},1);            % number of variational updates
Nx    = length(Ni);                     % number of variational updates
% Ne    = find([MDP.o(3,:) > 5,1],1) - 1; % how many outcomes were observed before scene choice (e.g. find how many saccades until first choice)
Ne = find([MDP.o(3,:) > 5,1],1);        % number of eye-movements until categorization/scene choice (including scene choice itself)

outcomes = {'null','UP','RIGHT','DOWN','LEFT'};
category_choices = {'UP RIGHT', 'RIGHT DOWN','DOWN LEFT','LEFT UP'};

    

for k = 1:Ne
    X(k,:) = x(MDP.o(3,k),:);
end

for j = 1:2
    T(:,j) = interp(X(:,j),8,2);
    T(:,j) = T(:,j) + spm_conv(randn(size(T(:,j))),2)/16;
end
       

for k = 1:Ne
    
    hold on;
    
    for i = 1:Nx
             
        if k < Ne
        
            for j = 1:4
                S{j} = zeros(size(UP));
            end
            
            for config_i = 1:size(all_configs,1)
                
                a = all_configs(config_i,:);
                a(strcmp(a,'n')) = outcomes(1);
                a(strcmp(a,'U')) = outcomes(2);
                a(strcmp(a,'R')) = outcomes(3);
                a(strcmp(a,'D')) = outcomes(4);
                a(strcmp(a,'L')) = outcomes(5);
                
                % posterior beliefs about scene
                %--------------------------------------------------
                p     = MDP.xn{1}(Ni(i),config_i,k,k);
                
                for j = 1:4
                    S{j} = S{j} + p.*eval(a{j});
                end
                
            end
            
            % image
            %------------------------------------------------------------------
            if i > 1
                delete(h);
            end
            for j = 1:numel(S)
                h(j) = imagesc(r + x(j + 1,1),r + x(j + 1,2),S{j}/max(S{j}(:)));
            end
            
%             stimulus
%             ------------------------------------------------------------------
            d   = (1 - exp(1 - i));
            
            switch MDP.o(1,k)
                case 1
                    imagesc(r,r,null.*mask*d)
                case 2
                    imagesc(r,r,UP.*mask*d)
                case 3
                    imagesc(r,r,RIGHT.*mask*d)
                case 4
                    imagesc(r,r,DOWN.*mask*d)
                case 5
                    imagesc(r,r,LEFT.*mask*d)
            end
            
        elseif k == Ne
            
            
            for j = 1:4
                S{j} = zeros(100,100,3);
            end
            
            for j = 1:4
                           
                % posterior beliefs about policy
                %--------------------------------------------------
                p     = MDP.un(j+5,k*Ni(i));
                
                S{j}(:,:,2:3) = S{j}(:,:,2:3) + p;
                
            end
            
            % image
            %------------------------------------------------------------------

            for j = 1:numel(S)
                imagesc(r + x(j + 5,1),r + x(j + 5,2),S{j});
            end
            
            
            % stimulus
            %------------------------------------------------------------------
            d   = (1 - exp(1 - i));
        
        
            switch MDP.o(2,k)
                case 1
                    imagesc(r,r,null.*mask*d)
                case 2
                    text(r,r,'CORRECT!')
                case 3
                    text(r,r,'INCORRECT!')
            end
        end
        
        for cat_i = 1:length(category_choices)
            text(x(cat_i+5,1),x(cat_i+5,2),category_choices{cat_i},'FontSize',14);
            hold on;
        end
        
        % save
        %------------------------------------------------------------------
        axis image ij; axis([-2,2.8,-2,2.8]),set(gca,'XColor','w','YColor','w'), drawnow
        M((k - 1)*Nx + i) = getframe(gca);
        
    end
    
    plot(X(k,1),X(k,2),'b.','MarkerSize',8);
    plot(T(1:(8*k),1),T(1:(8*k),1),'b:'); hold on;
    
    axis image; axis([-2,2.8,-2,2.8]), drawnow
    
        
end
    
%     % static pictures
%     %----------------------------------------------------------------------
%     subplot(2,Ne,Ne + k),hold off
%     for j = 1:numel(S)
%         imagesc(r + x(j + 1,1),r + x(j + 1,2),S{j}/max(S{j}(:))), hold on
%     end
%     
%     % stimulus
%     %------------------------------------------------------------------
%     switch MDP.o(1,k)
%         case 1
%             imagesc(r,r,null.*mask*d)
%         case 2
%             imagesc(r,r,UP.*mask*d)
%         case 3
%             imagesc(r,r,RIGHT.*mask*d)
%         case 4
%             imagesc(r,r,DOWN.*mask*d)
%         case 5
%             imagesc(r,r,LEFT.*mask*d)
%     end
%  
%     for j = 1:k
%         X(j,:) = x(MDP.o(3,j),:);
%     end
%     plot(X(:,1),X(:,2),'b.','MarkerSize',8)
%     
%     % save
%     %------------------------------------------------------------------
%     axis image, axis([-2,2,-2,2]), drawnow
    

% subplot(2,1,1)
% set(gca,'Userdata',{M,16})
% set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
% title('Scene construction','FontSize',16)
% title('Percept (click axis for movie)')
            
           
end

