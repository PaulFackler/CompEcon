% DEMDP03 Timber Cutting Model
  fprintf('\nDEMDP03 TIMBER CUTTING PROBLEM\n')
  close all  
  
% ENTER MODEL PARAMETERS
  price =  1;                                           % output price
  k     =  0.2;                                         % replanting cost
  sbar  =  0.5;                                         % carrying capacity
  gamma =  0.1;                                         % growth rate
  delta =  0.95;                                        % discount factor
  
% CONSTRUCT ACTION SPACE
  x = [0;1];            

% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfdp03';                                % reward/transition file
  model.discount = delta;                               % discount factor
  model.actions = x;                                    % model actions 
  model.params = {price k sbar gamma};                  % other parameters

% DEFINE APPROXIMATION SPACE
  n      = 500;                                         % degree of approximation
  fspace = fundefn('spli',n,0,sbar);                    % function space
  snodes = funnode(fspace);                             % state collocaton nodes

% SOLVE BELLMAN EQUATION  
  [c,s,v,x,resid] = dpsolve(model,fspace,snodes);
    
% OPTIMAL CUTTING STOCK LEVEL
  sstar=s(sum(x==0)+1);
  fprintf('Optimal Cutting Stock Level = %5.2f\n',sstar) 
  
% PLOT VALUE FUNCTION
  figure(1);
  plot(s,v);
  title('Value Function');
  xlabel('Biomass');
  ylabel('Value of Stand');
 
% PLOT SHADOW PRICE FUNCTION
  figure(2);
  p = funeval(c,fspace,s,1);
  plot(s,p);
  title('Shadow Price Function');
  xlabel('Biomass');
  ylabel('Marginal Value');

% PLOT RESIDUAL
  figure(3);
  plot(s,resid);
  title('Approximation Residual');
  xlabel('Biomass');
  ylabel('Residual');
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Timber Cutting Model')
