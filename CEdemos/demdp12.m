% DEMDP12 Production-Adjustment Model
  disp('DEMDP12 PRODUCTION ADJUSTMENT MODEL')
  close all 
  
% ENTER MODEL PARAMETERS
  alpha = 0.5;                              % marginal adjustment cost
  beta  = 0.5;                              % absolute demand elasticity
  kappa = 0.5;                              % unit cost of production
  sigma = 0.4;                              % demand shock standard deviation
  delta = 0.9;                              % discount factor
  
% COMPUTE SHOCK DISTRIBUTION
  m = 3;                                    % number of shocks
  [e,w] = qnwlogn(m,-sigma^2/2,sigma^2);             % shocks and proabilities

% PACK MODEL STRUCTURE  
  clear model
  model.func = 'mfdp12';                    % model functions
  model.discount = delta;                   % discount factor
  model.e = e;                              % shocks
  model.w = w;                              % probabilities
  model.params = {alpha beta kappa};        % other parameters
  
% COMPUTE CERTAINTY-EQUIVALENT STEADY-STATE
  estar = 1;                                % certainty equivalent shock
  xstar = ((1-beta)/kappa)^(1/beta);        % steady-state action
  sstar = [xstar 1];                        % steady-state states
  pstar = [xstar^(1-beta) alpha*(xstar-1)]; % steady-state shadow prices
  
% DEFINE APPROXIMATION SPACE
  n = [15 10];					                    % number of state collocation coordinates
  smin = [e(1) xstar-1.0];        					% minimum states
  smax = [e(m) xstar+3.0];							    % maximum states
  fspace = fundefn('cheb',n,smin,smax);     % function space
  scoord = funnode(fspace);	                % state collocation grid coordinates
  snodes = gridmake(scoord);			          % state collocation grid points
  
% CHECK MODEL DERIVATIVES AT C-E STEADY STATE
  dpcheck(model,sstar,xstar,estar);

% COMPUTE L-Q APPROXIMATION
  [vlq,xlq,plq,ss,xx,pp] = lqapprox(model,snodes,sstar,xstar,pstar);  

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',4);
  [c,s,v,x,resid] = dpsolve(model,fspace,snodes,vlq,xlq);
  
% COMPUTE SHADOW PRICES
  n = [length(s{1}) length(s{2})]; 
  p1 = funeval(c,fspace,s,[1 0]);  
  p2 = funeval(c,fspace,s,[0 1]);
  p1 = reshape(p1,n);     
  p2 = reshape(p2,n);
  
% PLOT OPTIMAL POLICY (Surface)
  figure(1)
  hh=surf(s{1},s{2},x');	
  title('Optimal Production Policy');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Production');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  
% PLOT OPTIMAL POLICY (Contours)
  figure(2)
  contour(s{1},s{2},x',25);
  title('Optimal Production Policy');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Production');
  rotate3d
  
% PLOT VALUE FUNCTION (Surface)
  figure(3)  
  hh=surf(s{1},s{2},v');
  title('Value Function');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
    
% PLOT VALUE FUNCTION (Contours)
  figure(4)
  contour(s{1},s{2},v',25);
  title('Value Function');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Optimal Value')
  rotate3d
    
% PLOT SHADOW PRICE FUNCTION 1 (Surface)
  figure(5)  
  hh=surf(s{1},s{2},p1');
  title('Shadow Price of Market Price');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Shadow Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  
% PLOT SHADOW PRICE FUNCTION 2 (Surface)
  figure(6)
  hh=surf(s{1},s{2},p2');
  title('Shadow Price of Lagged Production');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Shadow Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')

% PLOT SHADOW PRICE FUNCTION 1 (Contour)
  figure(7)
  contour(s{1},s{2},p1',25);  
  title('Shadow Price of Market Price');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Shadow Price')
  rotate3d
  
% PLOT SHADOW PRICE FUNCTION 2 (Contour)
  figure(8)
  contour(s{1},s{2},p2',25);
  title('Shadow Price of Lagged Production');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Shadow Price')
  rotate3d
  
% PLOT RESIDUAL
  figure(9)  
  hh=surf(s{1},s{2},resid');
  title('Approximation Residual');
  xlabel('Demand Shock');
  ylabel('Lagged Production'); 
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  
% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 20; nrep = 5000;
  sinit = smin(ones(nrep,1),:);
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);
  
% PLOT EXPECTED POLICY PATH
  figure(10);
  plot(0:nyrs,mean(xpath));
  title('Expected Policy Path');
  xlabel('Year');
  ylabel('Production');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Production-Adjustment Model',[1 3 9 10])