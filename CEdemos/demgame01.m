% DEMGAME01 Capital Production Game
  disp('DEMGAME01 CAPITAL PRODUCTION GAME')
  close all 
  
% ENTER MODEL PARAMETERS
  alpha = [8.0 4.0];
  beta  = [1.8 0.2];
  gamma = [0.4 3.0];
  psi   = 0.1;
  delta = 0.9;
  
% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfgame01';                  % model functions
  model.discount = delta;                   % discount factor
  model.params = {alpha,beta,gamma,psi};    % other parameters
 
% DEFINE APPROXIMATION SPACE
  n      = [8 8];                           % degree of approximation
  smin   = [0.7 0.7];                       % minimum state
  smax   = [1.3 1.3];                       % maximum state
  fspace = fundefn('cheb',n,smin,smax);     % function space
  scoord = funnode(fspace);                 % state collocaton coordinates
  snodes = gridmake(scoord);                % state collocaton nodes

% INITIALIZE POLICY, VALUE FUNCTIONS
  xinit = zeros(size(snodes));
  vinit = zeros(size(snodes));

  gamecheck(model,(smax+smin)/2,zeros(1,2));
  
% SOLVE BELLMAN EQUATIONS
  optset('gamesolve','nres',4);
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
  title('Investment: Player 1');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Investment');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([.7 1.3]);ylim([.7 1.3]); zlim([.04 .16]);
  view(-20,30);
  
% PLOT OPTIMAL POLICY
  figure(2)
  hh=surf(s{1},s{2},x(:,:,2)');  
  title('Investment: Player 2');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Investment');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([.7 1.3]);ylim([.7 1.3]);
  view(-80,30);
  
% PLOT VALUE FUNCTION
  figure(3)  
  hh=surf(s{1},s{2},v(:,:,1)');
  title('Value Function: Player 1');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([.7 1.3]);ylim([.7 1.3]);
  view(-20,30);
  
% PLOT VALUE FUNCTION
  figure(4)  
  hh=surf(s{1},s{2},v(:,:,2)');
  title('Value Function: Player 2');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Value');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([.7 1.3]);ylim([.7 1.3]);
  view(-80,30);
  
% PLOT OWN SHADOW PRICE FUNCTION
  figure(5)  
  hh=surf(s{1},s{2},p1');
  title('Own Shadow Price: Player 1');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([.7 1.3]);ylim([.7 1.3]);
  view(-20,30);
  
% PLOT OWN SHADOW PRICE FUNCTION
  figure(6)  
  hh=surf(s{1},s{2},p2');
  title('Own Shadow Price: Player 2');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Price')
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([.7 1.3]);ylim([.7 1.3]);
  view(-80,30);

% PLOT RESIDUAL
  figure(7)  
  hh=surf(s{1},s{2},resid(:,:,1)');
  title('Approximation Residual: Player 1');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([.7 1.3]);ylim([.7 1.3]); 
  view(-20,30);
  
% PLOT RESIDUAL
  figure(8)  
  hh=surf(s{1},s{2},resid(:,:,2)');
  title('Approximation Residual: Player 2');
  xlabel('S_1'); ylabel('S_2');
  zlabel('Residual');
  set(hh,'FaceColor','interp','EdgeColor','interp')
  xlim([.7 1.3]);ylim([.7 1.3]);
  view(-80,30);
  
% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 30;
  nrep =  2;
  sinit = smin(ones(nrep,1),:);
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);
  s1path = squeeze(spath(:,1,:));
  s2path = squeeze(spath(:,2,:));
  x1path = squeeze(xpath(:,1,:));

% PLOT EXPECTED STATE PATH
  figure(9);
  plot(0:nyrs,mean(s1path));
  title('Expected Capital Stock: Player 1');
  xlabel('Year');
  ylabel('Capital');

% PLOT EXPECTED STATE PATH
  figure(10);
  plot(0:nyrs,mean(s2path));
  title('Expected Capital Stock: Player 2');
  xlabel('Year');
  ylabel('Capital');

% PLOT EXPECTED POLICY PATH
  figure(11);
  plot(0:nyrs,mean(x1path));
  title('Expected Investment: Player 1');
  xlabel('Year');
  ylabel('Investment');

% PLOT EXPECTED POLICY PATH
  figure(12);
  plot(0:nyrs,mean(x1path));
  title('Expected Investment: Player 2');
  xlabel('Year');
  ylabel('Investment');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Capital-Production Game',[1 7 9 11])
