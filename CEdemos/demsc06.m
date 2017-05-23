% DEMSC06 Nonrenewable Resource Model (continuous time)
  close all 
  optset('scsolve','defaults'); 
  
  disp('DEMSC06 Nonrenewable Resource Model (continuous time)')
  
% Define Problem Parameters
  A=1;
  alpha=0.15;
  rho=0.05;

% Pack model structure
  clear model
  model.func='mfsc06';
  model.params={A,alpha,rho};
  
% Define nodes and basis matrices
  n=3;
  smin=0.01; smax=100;
  fspace=fundefn('cheb',n,smin,smax);
  snodes=funnode(fspace);

% Define initial values - identity mapping
  v0=snodes;
    
% Call solver 
  [cv,s,v,x,resid] = scsolve(model,fspace,snodes,v0);
 
% Plot optimal control function
  figure(1)
  plot(s,x)
  title('Optimal Extraction Rule')
  xlabel('S')
  ylabel('x')

% Plot value function
  figure(2)
  plot(s,v)
  title('Value Function')
  xlabel('S')
  ylabel('V(S)')

% Plot marginal value function
  figure(3)
  plot(s,funeval(cv,fspace,s,1))
  title('Marginal Value Function')
  xlabel('S')
  ylabel('V''(S)')

% Plot residual function
  figure(4)
  plot(s,resid)
  title('Approximation Residual')
  xlabel('S')
  ylabel('Residual')