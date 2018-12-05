function MDP_Drift_Diffusion_4outcomes

precis_colors = cool(length([1:0.2:2]));

color_iter = 1;

for i = 1:0.2:2                     % Change precision on A with each loop
    D{1} = [1 1 1 1]';
    
    %% interesting way to initialize likelihood matrix
%     a = repmat(spm_softmax(exp(i)*log([1,0,0,0]' + 4)),1,4);
%     shiftindex = [0:3];
%     [m,n]=size(a);
%     S=full(sparse(mod(shiftindex,m)+1,1:n,1,m,n));
%     a_new=ifft(fft(a).*fft(S),'symmetric');
        
    A{1} = spm_softmax(exp(i)*log(eye(4)+4));
    B{1}(:,:,1) = eye(4);           % Actions do nothing (the two B matrices here are only included to prevent the scheme entering 'HMM mode')
    B{1}(:,:,2) = eye(4);
    C{1} = zeros(4,1);
   
    MDP.A = A;
    MDP.B = B;
    MDP.C = C;
    MDP.D = D;
    MDP.T = 40;
    MDP.s = ones(1,MDP.T);
   
    mdp = spm_MDP_check(MDP);
    MDP = spm_MDP_VB_X(mdp);
    
    x = zeros(4,MDP.T);
    evidence_vecs = -ones(4) + 2*eye(4);
    for s_i = 1:4
        x(s_i,1) = evidence_vecs(s_i,:)*log(MDP.A{1}(MDP.o(1),:)');
    end
    
    for s_i = 1:4
        for t = 2:MDP.T
            x(s_i,t) = x(s_i,t-1) + evidence_vecs(s_i,:) * log(MDP.A{1}(MDP.o(t),:)');
        end
    end
%     x(1) = [1,-1]*log(MDP.A{1}(MDP.o(1),:)');
%     for t = 2:MDP.T
%         x(t) = x(t-1) + [1,-1]*log(MDP.A{1}(MDP.o(t),:)');
%     end

    for s_i = 1:4
        
        if s_i == MDP.s(1)
            plot(x(s_i,:),'Color',precis_colors(color_iter,:),'LineWidth',2,'DisplayName',sprintf('Evidence for true state, precision: %.1f',i));
        else
            plot(x(s_i,:),'Color',precis_colors(color_iter,:),'LineWidth',0.5,'DisplayName',sprintf('Evidence for other states, precision: %.1f',i));
        end
            
        hold on;
    end
    xlim([1 MDP.T])

%     plot(x), hold on
    clear MDP
    clear x
    color_iter = color_iter + 1;
end

legend('show')
xlabel('Sample')
ylabel('Evidence')
