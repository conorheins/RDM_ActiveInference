% prior beliefs about initial states (in terms of counts_: D and d
%--------------------------------------------------------------------------
d{1} = [1 1 1 1]';           % what:     {'A','B','C','D'}
d{2} = [1 0 0 0 0 0 0 0 0]'; % where:    {'start','1',...,'4','A','B','C',D'}
d{3} = [1 1]';               % flip(ud): {'no','yes'}
d{4} = [1 1]';               % flip(rl): {'no','yes'}
d{5} = [1 1]';               % flip (top right to lower right) : {'no','yes'}
d{6} = [1 1]';               % transpose (top right to lower left) : {'no','yes'}