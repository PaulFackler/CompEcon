%% DEMDDP08 Job Search Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
u     =  50;                            % weekly unemp. benefit
v     =  60;                            % weekly value of leisure
pfind = 0.90;                           % prob. of finding job
pfire = 0.10;                           % prob. of being fired
delta = 0.99;                           % discount factor

% State space
S = (1:2)';                             % vector of states
n = length(S);                          % number of states

% Action space (idle=1, active=2)
X = (1:2)';                             % vector of actions
m = length(X);                          % number of actions

% Reward function
f = zeros(n,m);
f(:,1) = v;                             % gets leisure
f(1,2) = u;                             % gets benefit

% State transition probability matrix
P = zeros(n,n,n);
P(1,:,1) = 1;                           % remains unemployed
P(2,1,1) = 1-pfind;                     % finds no job
P(2,1,2) = pfind;                       % finds job
P(2,2,1) = pfire;                       % gets fired
P(2,2,2) = 1-pfire;                     % keeps job

% Model structure
clear model
model.reward     = f;
model.transprob  = P;
model.discount   = delta;


%% SOLUTION

% Solve Bellman equation
wage   = (55:65)';
nw     = length(wage);
xtable = zeros(nw,2);
for i=1:nw
  f(2,2) = wage(i);
  model.reward = f;               % vary wage
  [v,x] = ddpsolve(model);        % solve via policy iteration
  xtable(i,:) = x';               % tabulate
end


%% ANALYSIS

% Display optimal policy
fprintf('\nOptimal Job Search Strategy')
fprintf('\n  (1=innactive, 2=active)\n')
fprintf('\nWage  Unemployed  Employed\n')
fprintf('%4i  %10i%10i\n',[wage xtable]')