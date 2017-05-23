% DEMDP04 American Option Pricing Model
  fprintf('\nDEMDP04 AMERICAN OPTION MODEL\n')
  close all  

% ENTER MODEL PARAMETERS
  sigma = 0.2;                                          % annual volatility
  T     = 0.5;                                          % years to expiration
  K     = 1.0;                                          % option strike price
  p0    = 1.0;                                          % current asset price 
  r     = 0.1;                                          % annual interest rate
  
% ENTER/COMPUTE TIME DISCRETIZATION PARAMETERS
  nt     = 30;
  deltat = T/(nt-1);
  beta   = exp(-r*deltat);

% COMPUTE SHOCK DISTRIBUTION
  m = 5;                                                % number of shocks
  [e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];
  nx = length(x);

% PACK MODEL STRUCTURE
  clear model
  model.horizon = nt;                                   % time horizon
  model.func = 'mfdp04';                                % model functions
  model.discount = beta;                                % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions   
  model.discretestates = 2;                             % index of discrete states 
  model.params={K};                                     % other parameters

% DEFINE APPROXIMATION SPACE
  d    = 6;                                             % state halfrange in std deviations
  ns   = 500;                                           % degree of approximation
  s0   = log(p0);                                       % current log of asset price
  smin = s0 - d*sigma*sqrt(T);                          % compute minimum state
  smax = s0 + d*sigma*sqrt(T);                          % compute maximum state
  fspace = fundefn('lin',ns,smin,smax,[],[0;1]);        % function space  
  scoord = funnode(fspace);                             % state collocaton nodes
  s = gridmake(scoord); 
  
% INITIALIZE VALUE FUNCTION
  p = exp(s(:,1));                                      % compute prices at states
  v = zeros(2*ns,1);                                    % terminal value function

% SOLVE BELLMAN EQUATION
  optset('dpsolve','nres',5);
  [c,s,v,x] = dpsolve(model,fspace,s,v);  

% PLOT VALUE FUNCTION
  figure(1);
  p = exp(scoord{1});                               
  v = reshape(v,ns,2,nt+1); 
  v = squeeze(v(:,1,:));
  plot(p,v(:,1),p,max(K-p,0),:);
  title('American Put Option Value');
  xlabel('Asset Price'); ylabel('Value');
 
% PLOT EXERCISE BOUNDARY
  figure(2);
  x = reshape(x,ns,2,nt);
  x = squeeze(x(:,1,:)); 
  plot(linspace(T,0,nt),p(sum(x==1)))
  title('American Put Option Optimal Exercise Boundary'); 
  xlabel('Time-to-Maturity'); ylabel('Asset Price');
    
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Soluton to the American Option Pricing Problem')