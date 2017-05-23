% DEMDP02 Industry Entry-Exit Model
  fprintf('\nDEMDP02 INDUSTRY ENTRY-EXIT MODEL\n')
  close all
  
% ENTER MODEL PARAMETERS
  pibar   = 0.0;                                        % long-run mean profits
  gamma  = 0.7;                                         % autoregressive coefficient
  kentry =  10;                                         % entry cost
  kexit  =   5;                                         % exit cost
  sigma  = 1.0;                                         % standard deviation of error
  delta  = 0.9;                                         % discount factor
  
% COMPUTE PRICE SHOCK DISTRIBUTION 
  m = 5;                                                % number of shocks
  [e,w] = qnwnorm(m,0,sigma.^2);
  
% CONSTRUCT ACTION SPACE
  x = [0;1];

% PACK MODEL STRUCTURE
  clear model  
  model.func = 'mfdp02';                                % model functions
  model.discount = delta;                               % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions
  model.discretestates = 2;                             % index of discrete states 
  model.params = {pibar gamma kentry kexit};            % other parameters

% DEFINE APPROXIMATION SPACE
  n = 500;                                              % number of profit nodes
  pimin = pibar+min(e)/(1-gamma)-.1;                    % minimum profit
  pimax = pibar+max(e)/(1-gamma)+.1;                    % maximum profit
  fspace = fundefn('spli',n,pimin,pimax,[],[0;1]);      % approximation space
  scoord = funnode(fspace);	                            % state collocation grid coordinates
  s = gridmake(scoord);	                                % state collocation grid points

% SOLVE BELLMAN EQUATION  
  [c,s,v,x,resid] = dpsolve(model,fspace,s);

% CRITICAL DECISION PRICES  
  sstar = s{1}(sum(x==0));
  fprintf('Critical Decision Prices\n')
  fprintf('   Enter   = %5.2f\n'  ,sstar(1)) 
  fprintf('   Exit    = %5.2f\n'  ,sstar(2))

% PLOT VALUE FUNCTION
  figure(1);
  plot(s{1},v);
  title('Value Function');
  legend('Inactive','Active',2);
  xlabel('Price');
  ylabel('Value');
  
% PLOT RESIDUAL
  figure(2);
  plot(s{1},resid);
  title('Approximation Residual');
  legend('Inactive','Active',2);
  xlabel('Price');
  ylabel('Residual');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename)
  prtfigs(mfilename,'Solution to the Industry Entry/Exit Model',[1 2])