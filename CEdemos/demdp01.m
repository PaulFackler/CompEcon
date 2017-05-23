% DEMDP01 Asset Replacement Model
function demdp01
  fprintf('\nDEMDP01 ASSET REPLACEMENT MODEL\n')
  close all
  
% ENTER MODEL PARAMETERS
  price   = 2.0;                                        % output price
  kbar    = 100;                                        % long-run mean replacement cost
  gamma   = 0.5;                                        % autoregressive coefficient
  abar    = 6;                                          % maximal age 
  sigma   = 15;                                         % standard deviation of price shock
  delta   = 0.9;                                        % discount factor  

% COMPUTE SHOCK DISTRIBUTION
  m = 5;                                                % number of shocks
  [e,w] = qnwnorm(m,0,sigma^2);                         % normal nodes and weights
  
% CONSTRUCT ACTION SPACE
  x = [0;1];

% PACK MODEL STRUCTURE
  clear model  
  model.func = 'mfdp01';                                % model functions
  model.discount = delta;                               % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.actions = x;                                    % model actions
  model.discretestates = 2;                             % index of discrete states 
  model.params = {price kbar gamma abar};               % other parameters

% DEFINE APPROXIMATION SPACE
  n = 100;                                              % number of cost nodes
  kmin =   20;                                          % minimum cost
  kmax =  200;                                          % maximum cost
  fspace = fundefn('spli',n,kmin,kmax,[],[1:abar]');    % approximation space
  snodes = funnode(fspace);	                            % state collocation grid coordinates
  s = gridmake(snodes);                                 % state collocation grid points


% CALL SOLVER
  [c,s,v,x,resid] = dpsolve(model,fspace,s);

% PREPARE FOR PLOTING
  n1 = length(s{1});
  n2 = length(s{2});
  v = reshape(v,n1,n2);
  x = reshape(x,n1,n2);
  resid = reshape(resid,n1,n2);
  
  sstar = s{1}(min(sum(x==1,1)+1,n1));
    
  legstr={'age=1'};
  for i=2:abar, legstr{i}=['age=' num2str(i)]; 
     if all(x(:,i)==1), break; end
  end
  ind=1:i;

% PRINT CRITICAL REPLACEMENT COST
  fprintf('\n Critical Replacement Cost \n');
  fprintf('     Age      Cost   \n');
  fprintf('   %5i    %6.0f   \n',[s{2}(ind) sstar(ind)]');

% PLOT VALUE FUNCTION
  figure(1)
  plot(s{1},v(:,ind));
  legend(legstr(ind),3);
  title('Value Function');
  xlabel('Replacement Cost');
  ylabel('Value');
  
% PLOT RESIDUAL
  figure(2);
  plot(s{1},resid(:,ind));
  legend(legstr(ind),1); 
  ylim([-.2 .35])
  title('Approximation Residual');
  xlabel('Replacement Cost');
  ylabel('Residual');
   
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to Asset Replacement Model',[1 2])
