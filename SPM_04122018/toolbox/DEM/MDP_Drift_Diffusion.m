function MDP_Drift_Diffusion
for i = 1:0.1:2                     % Change precision on A with each loop
    D{1} = [1 1]';
    A{1} = spm_softmax(exp(i)*log(eye(2)+4));
    B{1}(:,:,1) = eye(2);           % Actions do nothing (the two B matrices here are only included to prevent the scheme entering 'HMM mode')
    B{1}(:,:,2) = eye(2);
    C{1} = zeros(2,1);
   
    MDP.A = A;
    MDP.B = B;
    MDP.C = C;
    MDP.D = D;
    MDP.T = 40;
   
    mdp = spm_MDP_check(MDP);
    MDP = spm_MDP_VB_X(mdp);
    

    x(1) = [1,-1]*log(MDP.A{1}(MDP.o(1),:)');
    for t = 2:MDP.T
        x(t) = x(t-1) + [1,-1]*log(MDP.A{1}(MDP.o(t),:)');
    end
    
    
    plot(x), hold on
    clear MDP
    clear x
end