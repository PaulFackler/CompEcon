% DEMDDP06 Bioeconomic Model
  fprintf('\nDEMDDP06 BIOECONOMIC MODEL\n')
  close all

% Enter model parameters
  T     = 10;                   % foraging periods
  emax  =  8;                   % energy capacity
  e     = [2  4  4];            % energy from foraging
  p     = [1.0 0.7 0.8];        % predation survival probabilities
  q     = [0.5 0.8 0.7];        % foraging success probabilities

% Construct state space
  S = (0:emax)';                % energy levels
  n = length(S);                % number of states
  m = 3;                        % number of actions

% Construct reward function
  f = zeros(n,m);

% Construct state transition probability matrix
  P = zeros(m,n,n);
  for k=1:m
    P(k,1,1) = 1;
    for i=2:n
       % does not survive predation
       snext = 0;           j=getindex(snext,S); P(k,i,j) = P(k,i,j) + (1-p(k)); 
       % survives predation, finds food
       snext = S(i)-1+e(k); j=getindex(snext,S); P(k,i,j) = P(k,i,j) + p(k)*q(k);
       % survives predation, finds no food
       snext = S(i)-1;      j=getindex(snext,S); P(k,i,j) = P(k,i,j) + p(k)*(1-q(k));
    end
  end
   
% Initialize terminal value function
  vterm = ones(n,1);            % terminal value: survive
  vterm(1) = 0;                 % terminal value: death
  
% Pack model structure
  clear model
  model.reward     = f;
  model.transprob  = P;
  model.horizon    = T;
  model.discount   = 1;
  model.vterm      = vterm;

% Solve finite-horizon model using backward recursion
  [v,x] = ddpsolve(model);
  
% Plot survial probabilities, period 1
  figure(1); 
  h=bar(S,v(:,1),1);  set(h,'FaceColor',[.75 .75 .75])

  axis([-.5 emax+.5 0 1]);
  title('Survival Probability (Period 0)');
  xlabel('Stock of Energy'); ylabel('Probability');
  
% Plot survial probabilities, period 5
  figure(2); 
  h=bar(S,v(:,6),1);   set(h,'FaceColor',[.75 .75 .75])
  axis([-.5 emax+.5 0 1]);
  title('Survival Probability (Period 5)');
  xlabel('Stock of Energy'); ylabel('Probability');
  
% Save Plots as EPS Files
  prtfigs(mfilename,'Solution to the Bioeconomic Model',[1 2])