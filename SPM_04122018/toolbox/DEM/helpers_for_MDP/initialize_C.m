% priors: (utility) C
%--------------------------------------------------------------------------
T         = 6;
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{1}(6,:) =  2;                 % the agent expects to be right
C{1}(7,:) = -4;                 % and not wrong
C{2}(1:5,5:end) = -0.5;           % make tardy sampling costly , Berk said keep it 1:5 cuz 6,7,8,9 are the choice rows for the Scene identity
% C{2}(1:6,4:end) = -4;           % make tardy sampling costly 