% DEMREM01 Asset Pricing Model
  fprintf('\nDEMREM01 ASSET PRICING MODEL\n')
  close all  
  
% ENTER MODEL PARAMETERS
  dbar  = 1.0;                                          % long-run mean dividend
  gamma = 0.5;                                          % dividend autoregression coefficient
  beta  = 0.5;                                          % coefficient of risk aversion
  sigma = 0.1;                                          % dividend volatility
  delta = 0.9;                                          % discount factor

% COMPUTE SHOCK DISTRIBUTION
  m = 3;                                                % number of production shocks
  [e,w] = qnwnorm(m,0,sigma^2);                         % normal nodes and weights
  
% DEFINE APPROXIMATION SPACE
  n      = 10;                                          % degree of approximation
  dmin   = dbar+min(e)/(1-gamma);                       % minimum production
  dmax   = dbar+max(e)/(1-gamma);                       % maximum production
  fspace = fundefn('cheb',n,dmin,dmax);                 % function space
  dnode  = funnode(fspace);                             % state collocaton nodes

% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfrem01';                               % model functions
  model.discount = delta;                               % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.params = {delta dbar gamma beta};               % other parameters
  
% INITIALIZE RESPONSE
  xinit = dbar/(1-delta)+(dnode-dbar)/(1-gamma*delta); 
  
% SOLVE RATIONAL EXPECTATIONS EQULIBRIUM
  optset('remsolve','nres',10);
  optset('remsolve','showiters',0);
tic
  [c,d,p,x,f,resid] = remsolve(model,fspace,dnode,xinit); 
 toc
% PLOT EQUILIBRIUM PRICE
  figure(1); 
  plot(d,p);
  title('Equilibrium Pricing Function')
  xlabel('Dividend Level'); 
  ylabel('Asset Price');
  
% PLOT RESIDUAL
  figure(2); 
  plot(d,resid);
  title('Approximation Residual')
  xlabel('Dividend Level'); 
  ylabel('Residual');
  
% PLOT ARBITRAGE PROFIT
  figure(3); 
  plot(d,f);
  title('Arbitrage Profit Function')
  xlabel('Dividend Level'); 
  ylabel('Arbitrage Profit');
  
% SOLVE RATIONAL EXPECTATIONS EQULIBRIUM - DIRECT METHOD
  LHS = diag(dnode.^(-beta))*funbas(fspace,dnode);
  RHS = 0;
  for k=1:m
    dnext = dbar + gamma*(dnode-dbar) + e(k);
    LHS   = LHS - delta*w(k)*diag(dnext.^(-beta))*funbas(fspace,dnext);
    RHS   = RHS + delta*w(k)*dnext.^(1-beta);
  end
  c = LHS\RHS;

% COMPUTE RESIDUAL
  d = nodeunif(10*n,dmin,dmax);
  p = funeval(c,fspace,d);
  Eh=0;
  for k=1:m
    dnext = dbar + gamma*(d-dbar) + e(k);
    h     = diag(dnext.^(-beta))*(funeval(c,fspace,dnext)+dnext);
    Eh    = Eh + delta*w(k)*h;
  end
  resid = d.^(-beta).*funeval(c,fspace,d)-Eh;

% PLOT EQUILIBRIUM PRICE
  figure(4); 
  plot(d,p);
  title('Equilibrium Pricing Function')
  xlabel('Dividend Level'); 
  ylabel('Asset Price');
  
% PLOT RESIDUAL
  figure(5); 
  plot(d,resid);
  title('Approximation Residual')
  xlabel('Dividend Level'); 
  ylabel('Residual');

  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Asset Pricing Model',[4 5])
