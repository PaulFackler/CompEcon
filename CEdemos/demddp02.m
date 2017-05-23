% DEMDDP02 Asset Replacement Model
  fprintf('\nDEMDDP02 ASSET REPLACEMENT MODEL\n')
  close all
  
% Enter model parameters
  maxage  = 5;                          % maximum machine age
  repcost = 75;                         % replacement cost  
  delta   = 0.9;                        % discount factor  
 
% Construct state space
  S = (1:maxage)';                      % machine age
  n = length(S);                        % number of states

% Construct reward function (keep=1, replace=2)
  f = [50-2.5*S-2.5*S.^2 (50-repcost)*ones(n,1)]; 
  f(end,1) = -inf;
  
% Construct state transition function
  g = zeros(n,2);
  for i=1:n
     g(i,1) = min(i+1,n);                % keep
     g(i,2) = 1;                         % replace
  end
  
% Pack model structute
  clear model
  model.reward     = f;
  model.transfunc  = g;
  model.discount   = delta;

% Solve infinite-horizon model using policy iteration
  [v,x,pstar] = ddpsolve(model);

% Plot Optimal Value
  figure(1); plot(S,v); 
  title('Optimal Value Function');
  xlabel('Age of Machine'); ylabel('Value');
  
% Generate optimal path
  sinit = min(S); nyrs = 12;
  spath = ddpsimul(pstar,sinit,nyrs);
  
% Plot State Path
  figure(2); plot(0:nyrs,S(spath)); 
  title('Optimal State Path');
  xlabel('Year'); ylabel('Age of Machine');
  xlim([0 12])
  
% Save Plots as EPS Files
  prtfigs(mfilename,'Solution to the Asset Replacement Model',[1 2])