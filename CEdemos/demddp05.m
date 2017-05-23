% DEMDDP05 Water Management Model 
  fprintf('\nDEMDDP05 WATER MANAGEMENT MODEL\n')
  close all

% Enter model parameters
  alpha1 = 14;                      % producer benefit function parameter
  beta1  = 0.8;                     % producer benefit function parameter
  alpha2 = 10;                      % recreational user benefit function parameter
  beta2  = 0.4;                     % recreational user benefit function parameter  
  maxcap = 30;                      % maximum dam capacity
  r      = [0 1 2 3 4];             % rain levels
  p      = [0.1 0.2 0.4 0.2 0.1];   % rain probabilities
  delta  = 0.9;                     % discount factor
  
  evap=0;
  minlevel=10;

% Construct state space
  S = (0:maxcap)';              % vector of states
  n = length(S);                % number of states

% Construct action space
  X = (0:maxcap)';              % vector of actions
  m = length(X);                % number of actions

% Construct reward function
  f = zeros(n,m);
  for i=1:n
  for k=1:m
    if X(k)>S(i)
      f(i,k) = -inf;
    else
      f(i,k) = alpha1*X(k).^beta1+alpha2*max(0,S(i)-X(k)-minlevel).^beta2; 
      f(i,k) = alpha1*X(k).^beta1+alpha2*max(0,(S(i)-X(k)).*(S(i)-X(k)>=minlevel)).^beta2; 
    end
  end
  end
  
% Construct state transition matrix
  P = zeros(m,n,n);  
  for k=1:m
  for i=1:n
  for j=1:length(r)
    snext = max(0,min(S(i)-X(k)+r(j)-evap,maxcap));
    inext = getindex(snext,S);
    P(k,i,inext) = P(k,i,inext) + p(j);
  end
  end
  end

% Pack model structure
  clear model
  model.reward     = f;
  model.transprob  = P;
  model.discount   = delta;

% Solve infinite-horizon model using policy iteration
  [v,x,pstar] = ddpsolve(model);

% Plot optimal policy
  figure(1); 
  h=plot(S,X(x),'*'); % set(h,'FaceColor',[.75 .75 .75])
  axis([0 maxcap -inf inf]);
  title('Optimal Irrigation Policy');
  xlabel('Water Level'); ylabel('Irrigation');
  xlim([-1 31])
  ylim([0 6])

% Plot optimal value function
  figure(2); plot(S,v); 
  title('Optimal Value Function');
  xlabel('Water Level'); ylabel('Value');

% Generate random optimal paths, starting from zero water level
  sinit = ones(10000,1); 
  nyrs  = 30;
  spath = ddpsimul(pstar,sinit,nyrs);  

  % Plot expected water level over time, starting from zero water level
  figure(3)
  plot(0:nyrs,mean(S(spath)));
  title('Optimal State Path'); 
  xlabel('Year'); ylabel('Water Level');

% Compute steady state distribution of water level
  pi = markov(pstar);
  figure(4); 
  h=bar(S,pi,1); set(h,'FaceColor',[.75 .75 .75])
  title('Steady State Distribution');
  xlabel('Water Level'); ylabel('Probability'); 
  xlim([-1 31])
  ylim([0 0.16])

% Compute steady-state water level
  avgstock = pi'*S;
  fprintf('\nSteady-state Stock        %8.2f\n',avgstock)
  
% Save Plots as EPS Files
  prtfigs(mfilename,'Solution to the Water Management Model',[1 2 3 4])
