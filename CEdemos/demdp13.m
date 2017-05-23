% DEMDP13 Producton-Inventory Model
  fprintf('\nDEMDP13 PRODUCTION-INVENTORY MODEL\n')
  close all
  
% ENTER MODEL PARAMETERS
  c     = [0.5 0.1];                     % production cost function parameters
  k     = [0.1 0.1];                     % storage cost function parameters
  pbar  = 1.0;                           % long-run mean price
  rho   = 0.5;                           % mean-reversion coefficient
  sigma = 0.2;                           % standard deviation of price shocks
  delta = 0.9;                           % discount factor
  
% COMPUTE SHOCK DISTRIBUTION
  m   = 2;   					                   % number of shocks
  [e,w] = qnwnorm(m,0,sigma^2);          % shocks and proabilities
  
% PACK MODEL STRUCTURE  
  clear model
  model.func = 'mfdp13';                 % model functions
  model.discount = delta;                % discount factor
  model.e = e;                           % shocks
  model.w = w;                           % probabilities
  model.params = {c k pbar rho};         % other parameters
  
% COMPUTE CERTAINTY-EQUIVALENT STEADY-STATE   
  estar = 0;                             % certainty equivalent shock
  sstar = [0 pbar];                      % steady-state states    
  xstar = [(pbar-c(1))/c(2) 0];          % steady-state actions
  pstar = [pbar xstar(1)/(1-delta*rho)]; % steady-state shadow prices
  
% DEFINE APPROXIMATION SPACE
  n = [4 20];											       % number of state collocation coordinates
  smin = [0 pbar + min(e)/(1-rho)];			 % minimum states
  smax = [2 pbar + max(e)/(1-rho)];			 % maximum states
  fspace = fundefn('spli',n,smin,smax);	 % function space
  scoord = funnode(fspace);							 % state collocation grid coordinates
  snodes = gridmake(scoord);						 % state collocation grid points
  
% CHECK MODEL DERIVATIVES AT CE STEADY STATE
  dpcheck(model,sstar,xstar,estar);

% COMPUTE L-Q APPROXIMATION
  [vlq,xlq,plq,ss,xx,pp] = lqapprox(model,snodes,sstar,xstar,pstar);
  
% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  fprintf('\nHere we go!\n')
  [c,s,v,x,resid] = dpsolve(model,fspace,snodes,vlq,xlq);
  
% COMPUTE SHADOW PRICES
  n = [length(s{1}) length(s{2})]; 
  p1 = funeval(c,fspace,s,[1 0]);  
  p2 = funeval(c,fspace,s,[0 1]);
  p1 = reshape(p1,n);     
  p2 = reshape(p2,n);
  
% PLOT OPTIMAL POLICY (Surface)
  figure(1)
  hh=surf(s{1},s{2},x(:,:,1)');	
  title('Optimal Production Policy');
  xlabel('Beginning Inventories'); ylabel('Market Price');
  zlabel('Production');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  
% PLOT OPTIMAL POLICY (Surface)
  figure(2)
  hh=surf(s{1},s{2},x(:,:,2)');	
  title('Optimal Inventory Policy');
  xlabel('Beginning Inventories'); ylabel('Market Price');
  zlabel('Ending Inventories');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  
% PLOT VALUE FUNCTION (Surface)
  figure(3)  
  hh=surf(s{1},s{2},v');
  title('Value Function');
  xlabel('Beginning Inventories'); ylabel('Market Price');
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  
 % PLOT VALUE FUNCTION (Contours)
  figure(4)
  contour(s{1},s{2},v',25);
  title('Value Function');
  xlabel('Beginning Inventories'); ylabel('Market Price');
  zlabel('Optimal Value')
  rotate3d
  
% PLOT SHADOW PRICE FUNCTION 1 (Surface)
  figure(5)  
  hh=surf(s{1},s{2},p1');
  title('Shadow Price of Beginning Inventories');
  xlabel('Beginning Inventories'); ylabel('Market Price');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  
% PLOT SHADOW PRICE FUNCTION 1 (Contour)
  figure(6)
  contour(s{1},s{2},p1',25);
  title('Shadow Price of Beginning Inventories');
  xlabel('Beginning Inventories'); ylabel('Market Price');
  zlabel('Price')
  rotate3d
  
% PLOT RESIDUAL
  figure(7)  
  hh=surf(s{1},s{2},resid');
  title('Approximation Residual');
  xlabel('Beginning Inventories'); ylabel('Market Price');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  
% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 20;
  nrep = 10000;
  sinit = smin(ones(nrep,1),:);
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);
  s1path = squeeze(spath(:,1,:));
  s2path = squeeze(spath(:,2,:));
  x1path = squeeze(xpath(:,1,:));

% PLOT EXPECTED STATE PATH
  figure(8);
  plot(0:nyrs,mean(s1path));
  title('Expected State Path');
  xlabel('Year');
  ylabel('Beginning Inventories');

% PLOT EXPECTED STATE PATH
  figure(9);
  plot(0:nyrs,mean(s2path));
  title('Expected State Path');
  xlabel('Year');
  ylabel('Market Price');

% PLOT EXPECTED POLICY PATH
  figure(10);
  plot(0:nyrs,mean(x1path));
  title('Expected Policy Path');
  xlabel('Year');
  ylabel('Production');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Production-Inventory Model',[2 3 7 8])