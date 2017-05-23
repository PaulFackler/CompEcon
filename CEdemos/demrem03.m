% DEMREM03 Government Price Support Model
  fprintf('\nDEMREM03 GOVERNMENT PRICE SUPPORT MODEL\n')
  close all  

% ENTER MODEL PARAMETERS
  gamma = 2.0;                                          % inverse demand elasticity
  yvol  = 0.2;                                          % yield volatility
  xmax =  0.5;                                          % max storage
  c0 = 0.5;                                             % acreage supply, constant
  c1 = 0.5;                                             % acreage supply, slope   
  pstar = 0.9;                                          % government support price
  delta = 0.9;                                          % discount factor

% COMPUTE SHOCK DISTRIBUTION
  nshocks = 3;                                                % number of yield shocks
  [yshk,w] = qnwlogn(nshocks,0,yvol^2);                       % log normal nodes and weights
  
% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfrem03';                               % model functions
  model.e = yshk;                                       % shocks
  model.w = w;                                          % probabilities
  model.params = {delta xmax pstar gamma c0 c1};        % other parameters
    
% DEFINE APPROXIMATION SPACE
  n      = 100;                                          % degree of approximation
  smin   = 0.5;                                         % minimum supply
  smax   = 1.5;                                         % maximum supply
  fspace = fundefn('spli',n,smin,smax,[],yshk);         % function space
  scoord = funnode(fspace);	                            % state collocation grid coordinates
  snodes = gridmake(scoord);			                % state collocation grid points
     
% INITIALIZE STORAGE & PRICES
  nn = size(snodes,1);
  xinit = [zeros(nn,1) ones(nn,1)];                     % storage and acreage at nodes
  
% SOLVE RATIONAL EXPECTATIONS EQULIBRIUM
  optset('remsolve','nres',3);
  optset('arbit','lcpmethod','smooth');
  [c,s,x,p,f,resid] = remsolve(model,fspace,scoord,xinit); 

% PLOT GOVERNMENT STORAGE (Surface)
  figure(1)
  hh=surf(s{1},s{2},x(:,:,1)');	
  title('Government Policy');
  xlabel('Supply'); ylabel('Yield');
  zlabel('Government Storage');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  ylim([min(yshk) max(yshk)])
  zlim([0 0.45]);  
  
% PLOT PRIVATE ACREAGE (Surface)
  figure(2)
  hh=surf(s{1},s{2},x(:,:,2)');	
  title('Acreage Planted');
  xlabel('Supply'); ylabel('Yield');
  zlabel('Acreage');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  ylim([min(yshk) max(yshk)])
  zlim([0.7 1.1]);  
   
% PLOT PER ACRE REVENUE (Surface)
  figure(3)
  hh=surf(s{1},s{2},p');	
  title('Revenue');
  xlabel('Supply'); ylabel('Yield');
  zlabel('Revenue');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  ylim([min(yshk) max(yshk)])
    
% PLOT GOVERNMENT ARBITRAGE (Surface)
  figure(4)
  hh=surf(s{1},s{2},f(:,:,1)');	
  title('Government Arbitrage');
  xlabel('Supply'); ylabel('Yield');
  zlabel('Price Deficit');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  ylim([min(yshk) max(yshk)])
  
% PLOT PRIVATE ACREAGE ARBITRAGE (Surface)
  figure(5)
  hh=surf(s{1},s{2},f(:,:,2)');	
  title('Producer Arbitrage');
  xlabel('Supply'); ylabel('Yield');
  zlabel('Expected Profit');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  ylim([min(yshk) max(yshk)])
  
% PLOT RESIDUAL
  figure(6)  
  hh=surf(s{1},s{2},resid');
  title('Approximation Residual');
  xlabel('Supply'); ylabel('Yield');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  ylim([min(yshk) max(yshk)])
  zlim([-2 1.5]*1e-3);  
  
% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 10;
  nrep = 5000;
  sinit = [smax*ones(nrep,1) ones(nrep,1)];
  [spath,xpath] = remsimul(model,sinit,nyrs,s,x);
  s1path = squeeze(spath(:,1,:));
  x1path = squeeze(xpath(:,1,:));
  x2path = squeeze(xpath(:,2,:));

% PLOT EXPECTED STATE PATH
  figure(7);
  plot(0:nyrs,mean(s1path));
  title('Expected State Path');
  xlabel('Year');
  ylabel('Supply');
  
% PLOT EXPECTED POLICY PATH
  figure(8);
  plot(0:nyrs,mean(x1path));
  title('Expected Response Path');
  xlabel('Year');
  ylabel('Government Stocks');
    
% PLOT EXPECTED POLICY PATH
  figure(9);
  plot(0:nyrs,mean(x2path));
  title('Expected Response Path');
  xlabel('Year');
  ylabel('Acreage Planted');
  
% SAVE PLOTS AS EPS FILES
  prtfigs('demrem03','Solution to the Government Price Controls Model',[1 2 6 9])