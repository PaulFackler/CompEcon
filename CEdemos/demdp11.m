% DEMDP11 Monetary Policy Model
  disp('DEMDP11 MONETARY POLICY MODEL')
  close all 
  
% ENTER MODEL PARAMETERS
  alpha   = [0.9 0.4];              % transition function constant coefficients
  beta    = [0.8 0.5; 0.2 0.6];     % transition function state coefficients
  gamma   = [-0.8 0.0];             % transition function action coefficients
  omega   = [0.3 1];                % central banker's preference weights
  starget = [0 1];                  % equilibrium targets
  sigma   = 0.04*eye(2);            % shock covariance matrix 
  delta   = 0.9;                    % discount factor
  
% COMPUTE SHOCK DISTRIBUTION
  m   = [3 3];                    % number of shocks
  mu  = [0 0];                    % means of shocks
  [e,w] = qnwnorm(m,mu,sigma);    % shocks and proabilities

% PACK MODEL STRUCTURE  
  clear model
  model.func = 'mfdp11';                             % model functions
  model.discount = delta;                            % discount factor
  model.e = e;                                       % shocks
  model.w = w;                                       % probabilities
  model.params = {alpha beta gamma omega starget};   % other parameters
  
% DEFINE APPROXIMATION SPACE
  n = [10 10];					                 % number of state collocation coordinates
  nn = prod(n);                          % total number of collocation grid nodes
  smin = [-15 -10];                      % minimum states
  smax = [ 15  10];                      % maximum states
  fspace = fundefn('spli',n,smin,smax);  % function space
  scoord = funnode(fspace);	             % state collocation grid coordinates
  snodes = gridmake(scoord);		         % state collocation grid points
  
% COMPUTE CERTAINTY-EQUIVALENT STEADY-STATE   
  estar = mu;                                            % certainty equivalent shock
  sstar = starget;                                       % steady-state states    
  xstar = (sstar(1)-alpha(1)-beta(1,:)*sstar')/gamma(1); % steady-state actions
  pstar = [0 0];                                         % steady-state shadow prices
           
% CHECK MODEL DERIVATIVES AT CE STEADY STATE
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
  title('Optimal Monetary Policy');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Nominal Interest Rate');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([-15 15]); ylim([-10 10])
  
% PLOT OPTIMAL POLICY (Contours)
  figure(2)
  contour(s{1},s{2},x',25);
  title('Optimal Monetary Policy');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Nominal Interest Rate');
  rotate3d
  
% PLOT VALUE FUNCTION (Surface)
  figure(3)  
  hh=surf(s{1},s{2},v');
  title('Value Function');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([-15 15]); ylim([-10 10])
  
% PLOT VALUE FUNCTION (Contours)
  figure(4)
  contour(s{1},s{2},v',25);
  title('Value Function');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Optimal Value')
   
% PLOT SHADOW PRICE FUNCTION 1 (Surface)
  figure(5)  
  hh=surf(s{1},s{2},p1');
  title('Shadow Price of GDP Gap');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([-15 15]); ylim([-10 10])
  
% PLOT SHADOW PRICE FUNCTION 2 (Surface)
  figure(6)
  hh=surf(s{1},s{2},p2');
  title('Shadow Price of Inflation Rate');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([-15 15]); ylim([-10 10])

% PLOT SHADOW PRICE FUNCTION 1 (Contour)
  figure(7)
  contour(s{1},s{2},p1',25);
  title('Shadow Price of GDP Gap');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Price')
  
% PLOT SHADOW PRICE FUNCTION 2 (Contour)
  figure(8)
  contour(s{1},s{2},p2',25);
  title('Shadow Price of Inflation Rate');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Price')
  
% PLOT RESIDUAL
  figure(9)  
  hh=surf(s{1},s{2},resid');
  title('Approximation Residual');
  xlabel('GDP Gap'); 
  ylabel('Inflation Rate');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([-15 15]); ylim([-10 10]); zlim([-0.3 0.2]);
  

% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 20;
  nrep = 5000;
  sinit = smax(ones(nrep,1),:);
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);
  s1path = squeeze(spath(:,1,:));
  s2path = squeeze(spath(:,2,:));

% PLOT EXPECTED STATE PATH
  figure(10);
  plot(0:nyrs,mean(s1path));
  title('Expected State Path');
  xlabel('Year');
  ylabel('GDP Gap');

% PLOT EXPECTED STATE PATH
  figure(11);
  plot(0:nyrs,mean(s2path));
  title('Expected State Path');
  xlabel('Year');
  ylabel('Inflation Rate');

% PLOT EXPECTED POLICY PATH
  figure(12);
  plot(0:nyrs,mean(xpath));
  title('Expected Policy Path');
  xlabel('Year');
  ylabel('Interest Rate');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Monetary Policy Model',[1 9 10 11])