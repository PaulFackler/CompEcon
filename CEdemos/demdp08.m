% DEMDP08 Renewable Resource Model
  disp('DEMDP08 RENEWABLE RESOURCE MODEL')
  close all 
  
% ENTER MODEL PARAMETERS
  alpha = 4.0;                                          % growth function parameter
  beta  = 1.0;                                          % growth function parameter
  gamma = 0.5;                                          % demand function parameter
  kappa = 0.2;                                          % unit cost of harvest
  delta = 0.9;                                          % discount factor
    
% PACK MODEL STRUCTURE  
  clear model
  model.func = 'mfdp08';                                % model functions
  model.discount = delta;                               % discount factor
  model.params = {alpha beta gamma kappa};              % other parameters

% DEFINE APPROXIMATION SPACE
  n    = 8;                                             % degree of approximation
  smin = 6;                                             % minimum state
  smax = 9;                                             % maximum state
  fspace = fundefn('cheb',n,smin,smax);                 % function space
  snodes = funnode(fspace);                             % state collocaton nodes
  
% COMPUTE STEADY-STATE
  sstar = (alpha^2-1/delta^2)/(2*beta);                 % steady-state stock
  xstar = sstar - (delta*alpha-1)/(delta*beta);         % steady-state action
  pstar = xstar^(-gamma) - kappa;                       % steady-state shadow price    

% CHECK MODEL DERIVATIVES AT STEADY STATE
  dpcheck(model,sstar,xstar)

% COMPUTE L-Q APPROXIMATION
  [vlq,xlq] = lqapprox(model,snodes,sstar,xstar,pstar);  

% SOLVE BELLMAN EQUATION
  [c,s,v,x,resid] = dpsolve(model,fspace,snodes,vlq,xlq);
  
% COMPUTE L-Q APPROXIMATION FOR PLOTTING
  [vlq,xlq,plq] = lqapprox(model,s,sstar,xstar,pstar); 

% PLOT OPTIMAL POLICY
  figure(1);
  plot(s,x./s,s,xlq./s,sstar,xstar/sstar,'*');
  title('Optimal Harvest Policy');
  legend('Chebychev','L-Q');
  xlabel('Available Stock');
  ylabel('Harvest (% of Stock)');

% PLOT VALUE FUNCTION
  figure(2);
  plot(s,v,s,vlq)
  title('Value Function');
  legend('Chebychev','L-Q');
  xlabel('Available Stock');
  ylabel('Value');

% PLOT SHADOW PRICE FUNCTION
  figure(3);
  p = funeval(c,fspace,s,1);
  plot(s,p,s,plq,sstar,pstar,'*');
  title('Shadow Price Function');
  legend('Chebychev','L-Q');
  xlabel('Available Stock');
  ylabel('Price');

% PLOT RESIDUAL
  figure(4);
  plot(s,resid);
  title('Approximation Residual');
  xlabel('Available Stock');
  ylabel('Residual');
  ylim([-0.8 0.8]*1e-8)

% COMPUTE STATE AND POLICY PATH
  nyrs = 20;
  sinit = smin;
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);

% PLOT STATE PATH
  figure(5);
  plot(0:nyrs,spath);
  title('State Path');
  xlabel('Year');
  ylabel('Stock');

% PLOT POLICY PATH
  figure(6);
  plot(0:nyrs,xpath);
  title('Policy Path');
  xlabel('Year');
  ylabel('Harvest');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Renewable Resource Management Model',[1 3 4 5])