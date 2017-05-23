% DEMGAME02 Income Redistribution Game
  disp('DEMGAME02 INCOME REDISTRIBUTION GAME')
  close all  

% ENTER MODEL PARAMETERS
  alpha = [0.2 0.2];                        % utility parameter
  beta  = [0.5 0.5];                        % production elasticity
  gamma = [0.9 0.9];                        % capital survival rate
  psi   = 0.05;                             % risk sharing rate
  sigma = [0.1 0.1];                        % production shock volatility
  delta = 0.9;                              % discount factor

% COMPUTE SHOCK DISTRIBUTION
  nshocks = [3 3];                                % number of shocks
  [e,w] = qnwlogn(nshocks,[0 0],diag(sigma.^2));  % shocks and probabilities

% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfgame02';                  % model functions
  model.discount = delta;                   % discount factor
  model.e = e;                              % shocks
  model.w = w;                              % probabilities
  model.params = {alpha beta gamma psi};    % other parameters
    
% DEFINE APPROXIMATION SPACE
  n      = [15 15];                          % degree of approximation
  smin   = [ 3  3];                          % minimum state
  smax   = [11 11];                          % maximum state
  fspace = fundefn('cheb',n,smin,smax);      % function space
  scoord = funnode(fspace);                  % state collocaton coordinates
  snodes = gridmake(scoord);                 % state collocaton nodes

% INITIALIZE POLICY, VALUE, PRICE FUNCTIONS
  xstar=((1/delta-gamma)./beta).^(1./(beta-1));
  sstar=feval(model.func,'g',[0 0],xstar,w'*e,model.params{:});
  ratio=xstar./sstar;
  xinit = [snodes(:,1)*ratio(1) snodes(:,2)*ratio(2)];
  vinit=[];

  gamecheck(model,(smax+smin)/2,ratio.*(smax+smin)/2);
  
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
  title('Optimal Investment Policy Player 1');
  xlabel('Wealth 1'); ylabel('Wealth 2');
  zlabel('Investment');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([3 11]);ylim([3 11]);
  view(-45,30)
  
% PLOT OPTIMAL POLICY
  figure(2)
  hh=surf(s{1},s{2},x(:,:,2)');  
  title('Optimal Investment Policy Player 2');
  xlabel('Wealth 1'); ylabel('Wealth 2');
  zlabel('Investment');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([3 11]);ylim([3 11]);
  view(-45,30)
    
% PLOT VALUE FUNCTION
  figure(3)  
  hh=surf(s{1},s{2},v(:,:,1)');
  title('Value Function Player 1');
  xlabel('Wealth 1'); ylabel('Wealth 2');
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([3 11]);ylim([3 11]);
  view(-45,30)
  
% PLOT VALUE FUNCTION
  figure(4)  
  hh=surf(s{1},s{2},v(:,:,2)');
  title('Value Function Player 2');
  xlabel('Wealth 1'); ylabel('Wealth 2');
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([3 11]);ylim([3 11]);
  view(-45,30)
  
% PLOT OWN SHADOW PRICE FUNCTION
  figure(5)  
  hh=surf(s{1},s{2},p1');
  title('Own Shadow Price Player 1');
  xlabel('Wealth 1'); ylabel('Wealth 2');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([3 11]);ylim([3 11]);
  view(-45,30)
  
% PLOT OWN SHADOW PRICE FUNCTION
  figure(6)  
  hh=surf(s{1},s{2},p2');
  title('Own Shadow Price Player 2');
  xlabel('Wealth 1'); ylabel('Wealth 2');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([3 11]);ylim([3 11]);
  view(-45,30)

% PLOT RESIDUAL
  figure(7)  
  hh=surf(s{1},s{2},resid(:,:,1)');
  title('Approximation Residual Player 1');
  xlabel('Wealth 1'); ylabel('Wealth 2');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([3 11]);ylim([3 11]);
  view(-45,30)
  
% PLOT RESIDUAL
  figure(8)  
  hh=surf(s{1},s{2},resid(:,:,2)');
  title('Approximation Residual Player 2');
  xlabel('Wealth 1'); ylabel('Wealth 2');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([3 11]);ylim([3 11]);
  view(-45,30)

% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 20;
  nrep = 2000;
  sinit = 7.49*ones(nrep,2);
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);
  s1path = squeeze(spath(:,1,:));
  s2path = squeeze(spath(:,2,:));
  x1path = squeeze(xpath(:,1,:));
  x2path = squeeze(xpath(:,2,:));

% PLOT EXPECTED STATE PATH
  figure(9);
  plot(0:nyrs,mean(s1path));
  title('Expected Wealth Player 1');
  xlabel('Year');
  ylabel('Wealth');

% PLOT EXPECTED STATE PATH
  figure(10);
  plot(0:nyrs,mean(s2path));
  title('Expected Wealth Player 2');
  xlabel('Year');
  ylabel('Wealth');

% PLOT EXPECTED POLICY PATH
  figure(11);
  plot(0:nyrs,mean(x1path));
  title('Expected Investment Player 1');
  xlabel('Year');
  ylabel('Investment');

% PLOT EXPECTED POLICY PATH
  figure(12);
  plot(0:nyrs,mean(x2path));
  title('Expected Investment Player 2');
  xlabel('Year');
  ylabel('Investment');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Income Redistribution Game',[1 7 9 11])
