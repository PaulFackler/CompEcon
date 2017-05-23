% DEMREM02 Commodity Storage Model
  fprintf('\nDEMREM02 COMMODITY STORAGE MODEL\n')
  close all  

% ENTER MODEL PARAMETERS
  cost  = 0.1;                                          % unit storage cost
  gamma = 2.0;                                          % inverse demand elasticity
  yvol  = 0.2;                                          % yield volatility
  xmax =  1;                                          % max storage
  delta = 0.9;                                          % discount factor

% COMPUTE SHOCK DISTRIBUTION
  nshocks = 5;                                                % number of yield shocks
  [yshk,w] = qnwlogn(nshocks,0,yvol^2);                       % log normal nodes and weights
  
% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfrem02';                               % model functions
  model.e = yshk;                                       % shocks
  model.w = w;                                          % probabilities
  model.params = {delta,gamma,cost,xmax};               % other parameters
    
% DEFINE APPROXIMATION SPACE 
  n      = 50;                                          % degree of approximation
  smin   = min(yshk);                                   % minimum supply
  smax   = max(yshk)+xmax;                              % maximum supply
  fspace = fundefn('cheb',n,smin,smax);                 % function space
  snodes = funnode(fspace);                             % state collocaton nodes

% INITIALIZE STORAGE
  xinit = zeros(size(snodes));                          % storage at nodes
  
% SOLVE RATIONAL EXPECTATIONS EQULIBRIUM
  [c,s,x,h,f,resid] = remsolve(model,fspace,snodes,xinit);  

% PLOT EQUILIBRIUM STORAGE
  figure(1); 
  plot(s,x);
  title('Equilibrium Storage Function')
  xlabel('Supply'); 
  ylabel('Storage')
  ylim([0 0.9]);  


% PLOT EQUILIBRIUM PRICE
  figure(2); 
  plot(s,[h s.^-gamma]);
  title('Equilibrium Price Function')
  xlabel('Supply'); 
  ylabel('Price');
  ylim([0 3.5]);  
  
% PLOT ARBITRAGE PROFIT
  figure(3); 
  plot(s,f);
  title('Arbitrage Profit');
  xlabel('Supply');
  ylabel('Arbitrage Profit');

% PLOT RESIDUAL
  figure(4); 
  plot(s,resid);
  title('Approximation Residual');
  xlabel('Supply'); 
  ylabel('Residual');

% COMPUTE EXPECTED STOCK & STORAGE PATH
  nyrs = 20;
  nrep = 2000;
  sinit = smin*ones(nrep,1);
  [spath,xpath] = remsimul(model,sinit,nyrs,s,x);  

% PLOT EXPECTED STOCK PATH
  figure(5);
  plot(0:nyrs,mean(spath));
  title('Expected State Path');
  xlabel('Year'); 
  ylabel('Expected Supply');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Commodity Storage Model',[1 2 3 4])