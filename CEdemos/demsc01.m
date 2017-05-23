% DEMSC01 Optimal Growth Model
function demsc01
  close all  
  optset('scsolve','defaults'); 
  
  disp(' ')
  disp('DEMSC01 Optimal Growth Model')
  
% Define problem parameters
  alpha = 0.14;
  delta = 0.02;
  gamma = 0.5;
  rho   = 0.05;

% Pack model structure
  model.func='mfsc01';
  model.params={alpha,delta,gamma,rho};
  
% Define nodes and basis matrices
  n=20;
  smin=0.2; 
  smax=2;
  fspace=fundefn('cheb',n,smin,smax);
  snodes=funnode(fspace);

% Define initial values
  v0=((rho*snodes).^(1-gamma)-1)/(1-gamma)/rho;
  x0=(rho+delta)*snodes;
  
% Call solver  
  [cv,s,v,x,resid]=scsolve(model,fspace,snodes,v0);
  
% Get steady state
  [Kstar,Cstar]=ctsteadystate(model,fspace,cv);

% Display steady state results
  disp('Steady State Capital and Consumption')
  disp([Kstar Cstar])

% Plot value function
  figure(1)
  plot(s,v)
  title('Value Function')
  xlabel('K')
  ylabel('V(K)')
  xlim([0 2])

% Plot optimal control function
  figure(2)
  plot(s,x,Kstar,Cstar,'k*')
  title('Optimal Consumption Rule')
  xlabel('K')
  ylabel('C')
  xlim([0 2])

% Plot residual function
  figure(3)
  plot(s,resid)
  title('Approximation Residual')
  xlabel('K')
  ylabel('Residual')
  xlim([0 2])

% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Optimal Growth Model',[1 2 3])
