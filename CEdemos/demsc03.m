% DEMSC03 Production-Adjustment Model
  disp(' ')
  disp('DEMSC03 Production-Adjustment Model')
  close all 
  
% Define problem parameters
  kappa = 0.5;                                          % speed of price mean reversion
  alpha = 1;                                            % mean price
  sigma = 0.2;                                          % price volatility
  a     = 4;                                            % adjustment cost parameter
  c     = 2;                                            % variable cost parameter
  rho   = 0.1;                                          % discount rate  

% Pack model structure for computing long-run price density  
  clear model
  model.func = 'mfsc03b';                               % model functions
  model.params = {kappa alpha sigma a c rho};           % other parameters
  
% Evaluate long run price density 
  n = 50;                                          % number of state collocation coordinates
  pmin=0.4; pmax=2.5;                              % limits of approximation  
  fspace = fundefn('cheb',n,pmin,pmax);            % function space
  cp=itodensity(model,fspace);
  p=linspace(pmin,pmax,201)';
  figure(1)
  plot(p,funeval(cp,fspace,p));
  title('Long-run Price Density')
  xlabel('p')

% Pack model structure 
  clear model
  model.func = 'mfsc03';                        % model functions
  model.params = {kappa alpha sigma a c rho};   % other parameters
  
% Define nodes and basis matrices
  n = [10 10];					                        % number of state collocation coordinates
  smin = [pmin 0];                              % minimum states
  smax = [pmax 3];	                   		      % maximum states
  fspace = fundefn('cheb',n,smin,smax);         % function space
  scoord = funnode(fspace);	                    % state collocation grid coordinates
  
% Call solver 
  optset('scsolve','nres',3);
  [cv,s,v,x,resid] = scsolve(model,fspace,scoord);
  v=reshape(v,size(s{1},1),size(s{2},1));
  x=reshape(x,size(s{1},1),size(s{2},1));
  resid=reshape(resid,size(s{1},1),size(s{2},1));

% Compute shadow prices
  n = [length(s{1}) length(s{2})]; 
  p1 = funeval(cv,fspace,s,[1 0]);  
  p2 = funeval(cv,fspace,s,[0 1]);
  p1 = reshape(p1,n);     
  p2 = reshape(p2,n);
  
% Produce Plots  
% PLOT OPTIMAL POLICY (Surface)
  figure(2)
  hh=surf(s{1},s{2},x');	
  set(hh,'FaceColor','interp','EdgeColor','interp')
  title('Optimal Production Policy');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Production Adjustment Rate');
  view(-54,30)
  xlim([.4 2.5])
  
% PLOT OPTIMAL POLICY (Contours)
  figure(3)
  contour(s{1},s{2},x',25);
  title('Optimal Production Policy');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Production Adjustment Rate');
  
% PLOT VALUE FUNCTION (Surface)
  figure(4)  
  hh=surf(s{1},s{2},v');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  title('Value Function');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Value');
  view(-54,30)
  xlim([.4 2.5])
    
% PLOT VALUE FUNCTION (Contours)
  figure(5)
  contour(s{1},s{2},v',25);
  title('Value Function');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Value');
    
% PLOT SHADOW PRICE FUNCTION 1 (Surface)
  figure(6)  
  hh=surf(s{1},s{2},p1');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  title('Shadow Price of Market Price');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Shadow Price')
  view(-54,30)
  xlim([.4 2.5])
  
% PLOT SHADOW PRICE FUNCTION 2 (Surface)
  figure(7)
  hh=surf(s{1},s{2},p2');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  title('Shadow Price of Production Level');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Shadow Price')
  view(-54,30)

% PLOT SHADOW PRICE FUNCTION 1 (Contour)
  figure(8)
  contour(s{1},s{2},p1',25);
  title('Shadow Price of Market Price');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Shadow Price')
  xlim([.4 2.5])
  
% PLOT SHADOW PRICE FUNCTION 2 (Contour)
  figure(9)
  contour(s{1},s{2},p2',25);
  title('Shadow Price of Production Level');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Shadow Price')
  
% PLOT RESIDUAL
  figure(10)  
  hh=surf(s{1},s{2},resid');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  title('Approximation Residual');
  xlabel('Price'); 
  ylabel('Production Rate');
  zlabel('Residual');
  view(-54,30)
  xlim([.4 2.5])
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Production-Adjustment Model',[2 4 10 1])
