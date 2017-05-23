% DEMDDP08 Job Search Model
  fprintf('\nDEMDDP08 JOB SEARCH MODEL\n')
  close all

% Enter model parameters
  u     =  50;                  % weekly unemp. benefit
  v     =  60;                  % weekly value of leisure
  pfind = 0.90;                 % prob. of finding job
  pfire = 0.10;                 % prob. of being fired
  delta = 0.99;                 % discount factor

% Construct reward function
  f = zeros(2,2);
  f(:,1) = v;                   % gets leisure
  f(1,2) = u;                   % gets benefit

% Construct state transition probability matrix
  P        = zeros(2,2,2); 
  P(1,:,1) = 1;                 % remains unemployed
  P(2,1,1) = 1-pfind;           % finds no job
  P(2,1,2) = pfind;             % finds job
  P(2,2,1) = pfire;             % gets fired
  P(2,2,2) = 1-pfire;           % keeps job
  
% Pack model structure
  clear model
  model.reward     = f;
  model.transprob  = P;
  model.discount   = delta;

% Solve model at different wages
  xtable = [];
  wage   = 55:65;
  for w=wage
     f(2,2) = w; model.reward = f;  % vary wage
     [v,x]  = ddpsolve(model);      % solve via policy iteration
     xtable = [xtable x];           % tabulate
  end

% Print optimal policy
  fprintf('\nOptimal Job Search Strategy')
  fprintf('\n  (1=innactive, 2=active)\n')
  fprintf('\nWage  Unemployed  Employed\n')
  fprintf('%4i  %10i%10i\n',[wage;xtable])
  
% Save Plots as EPS Files
  prtfigs(mfilename)