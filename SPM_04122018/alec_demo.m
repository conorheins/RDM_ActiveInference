% the A matrix from the table I sent 
A = [1, 0, 0, 0;
     0, 1, 0, 0;
     0, 0, 1, 0;
     0, 0, 0, 1];

% make sure it is > 0 other the logarithm goes to Inf
A = A + exp(-16);
 
% take logarithm 
lnA = log(A);

% this is what we normally do
% agent recieves observation 1
o = 1;

% for simplicity ignore prior probability of hidden states
K = lnA(o,:);

% this is just the softmax function, ignore
K = (K - max(K));
Q = exp(K)/sum(exp(K));

% these are our posterior beliefs, as you can see, the agent correctly
% inferred it was looking at scene 1 location 1
disp(Q)

% now we want them to be uncertain whether they got observation 1 or 2,
% but certain it was location 1
o1 = 1;
o2 = 2;

% lets assume they distribution is 50% / 50%
K = lnA(o1,:) * 0.5 + lnA(o2,:) * 0.5;

% softmax
K = (K - max(K));
Q = exp(K)/sum(exp(K));

% as you can see, the agent knows it looking at location one, but is unsure
% whether it is scene 1 or scene 2. 
disp(Q)

% we can now look at its uncertainty about what it expects to observe at
% location 1
expected_o = A * Q';

% it is uncertain about whether it will recieve observation 1 or
% observation 2 at location 1, which makes sense
disp(expected_o)