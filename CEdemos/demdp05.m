% DEMDP05 Job Search Model
  fprintf('\nDEMDP05 JOB SEARCH MODEL\n')
  close all
  
% ENTER MODEL PARAMETERS
  u =  90;                                              % unemployment benefit
  v =  95;                                              % benefit of pure leisure
  pbar = 100;                                           % long-run mean wage
  gamma = 0.50;                                         % mean reversion rate
  pfind = 0.90;                                         % probability of finding job
  pkeep = 0.90;                                         % probability of keeping job
  delta = 0.95;                                         % discount factor
  sigma = 2;                                            % standard deviation of error
   
% COMPUTE SHOCK DISTRIBUTION 
  m = 7;                                                % number of shocks
  [e1,w1] = qnwnorm(m,0,sigma.^2);                      % price shocks
  e2 = [0;1]; w2 = [1-pfind;pfind];                     % job finding shock
  e3 = [0;1]; w3 = [1-pkeep;pkeep];                     % job keeping shock
  e = gridmake(e1,e2,e3);                               % shock vector
  w = ckron({w3,w2,w1});                                % shock probabilities

% CONSTRUCT ACTION SPACE
  x = [0;1];

% PACK MODEL STRUCTURE
  clear model  
  model.func = 'mfdp05';                                % model functions
  model.discount = delta;                               % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions
  model.discretestates = 2;                             % index of discrete states    
  model.params = {pbar gamma u v};                      % other parameters

% DEFINE APPROXIMATION SPACE
  n = 150;                                              % number of price nodes
  pmin = 60;                                            % minimum price
  pmax = 140;                                           % maximum price
  fspace = fundefn('spli',n,pmin,pmax,[],[0;1]);        % approximation space
  scoord = funnode(fspace);                             % state collocation grid coordinates
  s = gridmake(scoord);                                 % state collocation grid points

% SOLVE BELLMAN EQUATION  
  [c,s,v,x,resid] = dpsolve(model,fspace,s);

% PLOT VALUE FUNCTION
  figure(1);
  plot(s{1},v);
  title('Value Function');
  legend('Unemployed','Employed');
  xlabel('Wage Rate');
  ylabel('Value');
  
% PLOT RESIDUAL
  figure(2);
  plot(s{1},resid);
  title('Approximation Residual');
  legend('Unemployed','Employed');
  xlabel('Wage Rate');
  ylabel('Residual');
  
% CRITICAL DECISION PRICES  
  sstar = s{1}(sum(x==0));
  fprintf('Critical Decision Prices\n')
  fprintf('   Search if Unemployed = %5.2f\n'  ,sstar(1)) 
  fprintf('   Keep Job if Employed = %5.2f\n'  ,sstar(2))