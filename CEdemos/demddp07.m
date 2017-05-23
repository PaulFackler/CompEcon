% DEMDDP07 Renewable Resource Model
  fprintf('\nDEMDDP07 RENEWABLE RESOURCE MODEL\n')
  close all

% Enter model parameters
  delta =  0.9;                  % discount factor
  alpha =  4.0;                  % growth function parameter
  beta  =  1.0;                  % growth function parameter
  gamma =  0.5;                  % demand function parameter
  cost  =  0.2;                  % unit cost of harvest

% Construct state space
  n    = 200;                    % number of states
  smin = 0;                      % minimum state
  smax = 8;                      % maximum state
  S    = nodeunif(n,smin,smax);  % vector of states

% Construct action space
  m    = 100;                    % number of actions
  xmin = 0;                      % minimum action
  xmax = 6;                      % maximum action
  X    = nodeunif(m,xmin,xmax);  % vector of actions

% Construct reward function
  f = zeros(n,m);
  for k=1:m
     f(:,k) = (X(k).^(1-gamma))/(1-gamma)-cost*X(k);
     f(S<X(k),k) = -inf;
  end

% Construct state transition function
  g = zeros(n,m);
  for i=1:n
  for k=1:m
    snext  = alpha*(S(i)-X(k)) - 0.5*beta*(S(i)-X(k)).^2;
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
  title('Optimal Harvest Policy');
  xlabel('Stock'); ylabel('Harvest');

% Plot optimal value
  figure(2); plot(S,v); 
  title('Optimal Value Function');
  xlabel('Stock'); ylabel('Value');
    
% Generate optimal path starting at the maximal state level
  nyrs = 20;
  spath = ddpsimul(pstar,n,nyrs);
  
% Plot stock levels over time starting from carrying capacity stock level
  figure(3);
  plot(0:nyrs,S(spath)); 
  xlabel('Year'); ylabel('Stock');
  title('Optimal State Path')

% Plot optimal transition function
  figure(4)
  [ii,jj]=find(pstar);
  plot(S(ii),S(jj),S,S,'--');
  xlabel('S(t)'); ylabel('S(t+1)');
  title('Optimal State Transition Function')
  
% Save Plots as EPS Files
  prtfigs(mfilename)
