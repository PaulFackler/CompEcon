% DEMDP09 Non-Renewable Resource Model
  disp('DEMDP09 NON-RENEWABLE RESOURCE MODEL')
  close all  

% ENTER MODEL PARAMETERS
  a = [10 0.8];
  b = [12 1.0];
  delta = 0.9;

% PACK MODEL STRUCTURE
  model.func = 'mfdp09';
  model.discount = delta;
  model.params = {a b};

% DEFINE APPROXIMATION SPACE
  n      = 101;
  smin   =   0;
  smax   =  10;
  fspace = fundefn('spli',n,smin,smax);
  snodes = funnode(fspace);

% INITIALIZE POLICY, VALUE, PRICE  FUNCTIONS
  xinit = zeros(size(snodes));                  % initial policy function
  vinit = zeros(size(snodes));                  % initial value function

% CHECK MODEL DERIVATIVES AT INITIAL VALUES
  dpcheck(model,snodes,xinit);
  
% SOLVE BELLMAN EQUATION
  [c,s,v,x,resid] = dpsolve(model,fspace,snodes,vinit,xinit);

% PLOT OPTIMAL POLICY
  figure(1);
  plot(s,x,s,s,':');
  title('Optimal Harvest Policy');
  xlabel('Available Stock');
  ylabel('Harvest');

% PLOT VALUE FUNCTION
  figure(2);
  plot(s,v);
  title('Value Function');
  xlabel('Available Stock');
  ylabel('Value');

% PLOT SHADOW PRICE FUNCTION
  figure(3);
  p = funeval(c,fspace,s,1);
  plot(s,p);
  title('Shadow Price Function');
  xlabel('Available Stock');
  ylabel('Price');
  ylim([-1 7])

% PLOT RESIDUAL
  figure(4);
  plot(s,resid);
  title('Approximation Residual');
  xlabel('Available Stock');
  ylabel('Residual');

% COMPUTE STATE AND POLICY PATH
  nyrs = 20;
  sinit = smax;
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
  prtfigs(mfilename,'Solution to the Nonrenewable Resource Management Model',[1 3 4 5])