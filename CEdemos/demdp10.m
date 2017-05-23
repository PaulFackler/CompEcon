% DEMDP10 Water Management Model
  disp('DEMDP10 WATER MANAGEMENT MODEL')
  close all  

% ENTER MODEL PARAMETERS
  a = [ 1  2];                           % producer benefit function parameter
  b = [-2 -3];                           % recreational user benefit function parameter
  sigma = 0.2;                           % rainfall standard deviation
  delta = 0.9;                           % discount factor

% COMPUTE GAUSSIAN NODES AND WEIGHTS
  m = 3;                                 % number of shocks
  [e,w] = qnwlogn(m,-sigma^2/2,sigma^2);          % shocks and proabilities
  
  % DEFINE APPROXIMATION SPACE
  n    = 100;                             % degree of approximation
  smin =  2;                             % minimum state
  smax =  4.5;                             % maximum state
  fspace = fundefn('lin',n,smin,smax);  % function space
  snodes = funnode(fspace);              % state collocaton nodes

% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfdp10';                 % model function file
  model.discount = delta;                % discount factor
  model.e = e;                           % shocks
  model.w = w;                           % probabilities
  model.params = {a b smax};                  % other parameters
  
% COMPUTE CERTAINTY-EQUIVALENT STEADY-STATE
  estar = 1;                                    % certainty equivalent shock
  xstar = 1;                                    % steady-state action
  sstar = 1+(a(1)*(1-delta)/a(2))^(1/b(2));     % steady-state stock
  pstar = a(1);                                 % steady-state shadow price  

% CHECK MODEL DERIVATIVES AT CE STEADY STATE
  dpcheck(model,sstar,xstar,estar) 

% COMPUTE L-Q APPROXIMATION
  [vlq,xlq] = lqapprox(model,snodes,sstar,xstar,pstar);  
  
% SOLVE BELLMAN EQUATION
  [c,s,v,x,resid] = dpsolve(model,fspace,snodes,vlq,xlq);
  
% COMPUTE L-Q APPROXIMATION FOR PLOTTING 
  [vlq,xlq,plq] = lqapprox(model,s,sstar,xstar,pstar);
  
% PLOT OPTIMAL POLICY
  figure(1);
  plot(s,x,s,xlq,sstar,xstar,'*');
  title('Optimal Irrigation Policy');
  legend('Chebychev','L-Q',2);
  xlabel('Water Level');
  ylabel('Irrigation');

% PLOT VALUE FUNCTION
  figure(2);
  plot(s,v,s,vlq)
  title('Value Function');
  legend('Chebychev','L-Q');
  xlabel('Water Level');
  ylabel('Value');

% PLOT SHADOW PRICE FUNCTION
  figure(3);
  p = funeval(c,fspace,s,1);
  plot(s,p,s,plq,sstar,pstar,'*');
  title('Shadow Price Function');
  legend('Chebychev','L-Q');
  xlabel('Water Level');
  ylabel('Price');

% PLOT RESIDUAL
  figure(4);
  plot(s,resid);
  title('Approximation Residual');
  xlabel('Water Level');
  ylabel('Residual');

% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 30;
  nrep = 10000;
  sinit = smin*ones(nrep,1);
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);  

% PLOT EXPECTED STATE PATH
  figure(5);
  plot(0:nyrs,mean(spath));
  title('Expected State Path');
  xlabel('Year');
  ylabel('Water Level');
  
% PLOT EXPECTED POLICY PATH
  figure(6);
  plot(0:nyrs,mean(xpath));
  title('Expected Policy Path');
  xlabel('Year');
  ylabel('Irrigation');
  
% COMPUTE STEADY-STATE DISTRIBUTION
  nsmooth = 7; nbin = 80;
  [ss,pi,xx] = dpstst(model,nsmooth,nbin,s,x);
  sstar = pi'*ss;  
  xstar = pi'*xx;

% PLOT STEADY-STATE DISTRIBUTION
  figure(7);
  h=bar(ss,pi,1); set(h,'FaceColor',[.75 .75 .75])
  title('Steady State Distribution');
  xlabel('Water Level');
  ylabel('Probability');
  ylim([0 0.08])  

  fprintf('Steady State Means\n') 
  fprintf('   Stock        = %5.2f\n'  ,sstar)
  fprintf('   Irrigation   = %5.2f\n'  ,xstar)
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Water Management Model',[1 3 4 7])