% DEMGAME03 Marketing Board Game
  disp('DEMGAME03 MARKETING BOARD GAME')
  close all 
  
% ENTER MODEL PARAMETERS
  kappa = 0.05;                           % unit storage cost
  gamma = -0.5;                           % inverse demand elasticity
  xmax  = [0.2 0.2];                      % maximum storage
  mu    = log([0.5 0.5]);                 % relative size of agents
  yvol  = 0.2;                            % yield volatility
  delta = 0.95;                           % discount factor
  
% COMPUTE SHOCK DISTRIBUTION
  nshocks = [3 3];                        % number of shocks
  cov     = (yvol^2)*eye(2);              % covariance matrix of log shocks
  [e,w]   = qnwlogn(nshocks,mu,cov);      % shocks and proabilitities
  
% PACK MODEL STRUCTURE
  model.func = 'mfgame03';                % model functions
  model.discount = delta;                 % discount factor
  model.e = e;                            % shocks
  model.w = w;                            % probabilities
  model.params={kappa gamma xmax};        % other parameters

% DEFINE APPROXIMATION SPACE
  n      = [20 20];                       % degree of approximation
  smin   = min(e);                        % minimum supply
  smax   = max(e)+xmax;                   % maximum supply
  fspace = fundefn('spli',n,smin,smax);   % function space
  scoord = funnode(fspace);               % state collocaton coordinates
  snodes = gridmake(scoord);              % state collocaton nodes

% INITIALIZE POLICY, VALUE FUNCTIONS
  xinit = zeros(size(snodes));            % initial policy functions
  vinit = zeros(size(snodes));            % initialize value function

  gamecheck(model,(smax+smin)/2,zeros(1,2));
  
% SOLVE BELLMAN EQUATIONS
  optset('gamesolve','nres',3);
  [c,s,v,x,resid] = gamesolve(model,fspace,vinit,xinit);
  
% COMPUTE SHADOW PRICES
  n  = [length(s{1}) length(s{2})]; 
  p1 = funeval(c(:,1),fspace,s,[1 0]);  
  p2 = funeval(c(:,2),fspace,s,[0 1]);
  p1 = reshape(p1,n);     
  p2 = reshape(p2,n);
  
% PLOT OPTIMAL POLICY 
  figure(1)
  hh=surf(s{1},s{2},x(:,:,1)');  
  title('Optimal Storage: Country 1');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Storage');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])
  zlim([0 .12]);
  
% PLOT OPTIMAL POLICY
  figure(2)
  hh=surf(s{1},s{2},x(:,:,2)');  
  title('Optimal Storage: Country 2');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Storage');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])
  
% PLOT VALUE FUNCTION
  figure(3)  
  hh=surf(s{1},s{2},v(:,:,1)');
  title('Value Function: Country 1');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])
  
% PLOT VALUE FUNCTION
  figure(4)  
  hh=surf(s{1},s{2},v(:,:,2)');
  title('Value Function: Country 2');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])
  
% PLOT OWN SHADOW PRICE FUNCTION
  figure(5)  
  hh=surf(s{1},s{2},p1');
  title('Own Shadow Price: Country 1');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])
  
% PLOT OWN SHADOW PRICE FUNCTION
  figure(6)  
  hh=surf(s{1},s{2},p2');
  title('Own Shadow Price: Country 2');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])

% PLOT RESIDUAL
  figure(7)  
  hh=surf(s{1},s{2},resid(:,:,1)');
  title('Approximation Residual: Country 1');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])
  
% PLOT RESIDUAL
  figure(8)  
  hh=surf(s{1},s{2},resid(:,:,2)');
  title('Approximation Residual: Country 2');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])
  
% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 10;
  nrep = 1000;
  sinit = smax(ones(nrep,1),:);
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);
  s1path = squeeze(spath(:,1,:));
  s2path = squeeze(spath(:,2,:));
  x1path = squeeze(xpath(:,1,:));
  x2path = squeeze(xpath(:,2,:));

% PLOT EXPECTED STATE PATH
  figure(9);
  plot(0:nyrs,mean(s1path));
  title('Expected Supply: Country 1');
  xlabel('Year');
  ylabel('Supply');

% PLOT EXPECTED STATE PATH
  figure(10);
  plot(0:nyrs,mean(s2path));
  title('Expected Supply: Country 2');
  xlabel('Year');
  ylabel('Supply');

% PLOT EXPECTED POLICY PATH
  figure(11);
  plot(0:nyrs,mean(x1path));
  title('Expected Storage: Country 1');
  xlabel('Year');
  ylabel('Storage');

% PLOT EXPECTED POLICY PATH
  figure(12);
  plot(0:nyrs,mean(x2path));
  title('Expected Storage: Country 2');
  xlabel('Year');
  ylabel('Storage');
  
% PLOT MARKET PRICE FUNCTION
  figure(13)  
  ss = gridmake(s);
  xx = reshape(x,size(ss));
  pp = sum(ss-xx,2).^gamma;
  pp = reshape(pp,n); 
  hh=surf(s{1},s{2},pp');
  title('Market Price');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Market Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([s{1}(1),s{1}(end)])
  ylim([s{2}(1),s{2}(end)])
  zlim([0.7 1.2]);  
  view(-37,52)
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Marketing Board Game',[1 7 13 9])
