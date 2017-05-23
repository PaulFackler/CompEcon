% DEMDDP03 Asset Replacement with Maintenance Model
  fprintf('\nDEMDDP03 ASSERT REPLACEMENT WITH MAINTENANCE MODEL\n')
  close all

% Enter model parameters
  maxage  = 5;                          % maximum asset age
  repcost = 75;                         % replacement cost  
  mancost = 10;                         % maintenance cost
  delta   = 0.9;                        % discount factor  
 
% Construct state space
  s1 = (1:maxage)';                     % asset age
  s2 = (0:maxage-1)';                   % servicings
  S  = gridmake(s1,s2);                 % combined state grid
  n  = length(S);                       % total number of states

% Construct reward function (no action=1, service=2, replace=3)
  q  = (50-2.5*S(:,1)-2.5*S(:,1).^2);
  f1 = q.*min(1,1-(S(:,1)-S(:,2))/maxage);
  f2 = q.*min(1,1-(S(:,1)-S(:,2)-1)/maxage) - mancost;
  f3 = (50 - repcost)*ones(n,1);
  f  = [f1 f2 f3];
  
% Construct state transition function
  g = zeros(n,3);
  for i=1:n
    g(i,1) = getindex([S(i,1)+1 S(i,2)],S);
    g(i,2) = getindex([S(i,1)+1 S(i,2)+1],S);
    g(i,3) = getindex([1 0],S);
  end
  
% Pack model structure
  clear model
  model.reward     = f;
  model.transfunc  = g;
  model.discount   = delta;

% Solve infinite-horizon model using policy iteration
  [v,x,pstar] = ddpsolve(model);
  
% Generate optimal path
  sinit = 1; 
  nyrs  = 12;
  [spath,xpath] = ddpsimul(pstar,sinit,nyrs,x); 
  
% Plot State Path (Age)
  figure(1); plot(0:nyrs,S(spath,1)); 
  title('Optimal State Path');
  xlabel('Year'); ylabel('Age of Asset');
  xlim([0 12])

  
% Plot State Path (Servicings)
  figure(2); plot(0:nyrs,S(spath,2)); 
  title('Optimal State Path');
  xlabel('Year'); ylabel('Number of Servicings');
  xlim([0 12])
  ylim([0 2.25])
  
% Save Plots as EPS Files
  prtfigs(mfilename,'Solution to the Asset Replacement-Maintenance Model',[1 2])