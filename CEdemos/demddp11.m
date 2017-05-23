% DEMDDP11 Optimal growth model
  fprintf('\nDEMDDP11 OPTIMAL GROWTH MODEL\n')
  close all

% Enter model parameters
  delta =  0.9;                 % discount factor
  alpha =  0.2;                 % utility parameter
  beta  =  0.5;                 % production parameter
  gamma =  0.9;                 % capital survival rate
  sigma =  0.1;                 % production shock volatility

% Enter program parameters
  smin  =  1.0;                 % minimum state
  smax  = 10.0;                 % maximum state

% Construct state space
  n =  200;                     % number of states
  S = nodeunif(n,smin,smax);    % vector of states

% Construct action space
  m =  500;                     % number of actions
  X = nodeunif(m,0.65,0.99);    % vector of actions

% Construct reward function
  f = zeros(n,m);
  for k=1:m
     f(:,k) = ((S-X(k)*S).^(1-alpha))/(1-alpha);
  end

% Construct state transition matrix
  g = zeros(n,m);
  for i=1:n
  for k=1:m
    snext = gamma*X(k)*S(i) + (X(k)*S(i)).^beta;
    g(i,k) = getindex(snext,S);
  end
  end

% Pack model structure
  clear model
  model.reward     = f;
  model.transfunc  = g;
  model.discount   = delta;

% Solve infinite-horizon model using policy iteration
  [v,x,pstar] = ddpsolve(model);

% Plot optimal policy
  figure(1); plot(S,X(x));
  title('Optimal Investment');
  xlabel('Wealth'); ylabel('Investment as Precent of Wealth');

% Plot optimal value
  figure(2); plot(S,v);
  title('Optimal Value Function'); 
  xlabel('Wealth'); ylabel('Value');
  
% Generate optimal path
  nyrs = 20;
  st = ddpsimul(pstar,smin,nyrs);

% Plot stock levels over time starting from stocks = smax
  figure(3); plot(0:nyrs,S(st));
  title('Optimal State Path'); 
  xlabel('Year'); ylabel('Wealth');

% Compute steady state distribution
  pi = markov(pstar);

% Compute steady-state stock
  avgstock = pi'*S;
  fprintf('\nSteady-state Wealth     %8.2f\n',avgstock)
  
% Save Plots as EPS Files
  prtfigs(mfilename)
