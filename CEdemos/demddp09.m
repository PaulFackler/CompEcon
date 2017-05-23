% DEMDDP09 Deterministic Cow Replacement Model
  fprintf('\nDEMDDP09 DETERMINISTIC COW REPLACEMENT MODEL\n')
  close all

% Enter model parameters
  delta = 0.9;                  % discount factor  
  cost  = 500;                  % replacement cost  
  price = 150;                  % milk price

% Construct state space
  S = (1:10)';                  % lactation states
  n = length(S);                % number of states

% Construct action space
  X = ['K';'R'];                % keep or replace
  m = length(X);                % number of actions

% Construct reward function (actions keep=1, replace=2)
  y = (-0.2*S.^2+2*S+8);        % yield per lactation  
  f = [price*y price*y-cost];   % net revenue by action
  f(10,1) = -inf;               % force replace at lactation 10

% Construct state transition matrix
  g = zeros(n,m);
  for i=1:n
     g(i,1) = min(i+1,n);       % Raise lactation number by 1, if keep
     g(i,2) = 1;                % Lactation number reverts to 1, if replace
  end
  
% Pack model structure
  clear model
  model.reward     = f;
  model.transfunc  = g;
  model.discount   = delta;

% Solve infinite-horizon model using policy iteration
  [v,x] = ddpsolve(model);

% Plot optimal policy and value function
  figure(1); bar(S,x);  xlabel('Age'); ylabel('Optimal Decision');
  figure(2); plot(S,v); xlabel('Age'); ylabel('Optimal Value');

% Display optimal policy and value function
  disp('Lactation, optimal policy, value function')
  fprintf('%2i %2i %8.1f\n',[S x  v]')
  
% Save Plots as EPS Files
  prtfigs(mfilename)