% DEMDP06 Asset Replacement-Maintanence Model
  fprintf('\nDEMDP06 ASSET REPLACEMENT-MAINTANENCE MODEL\n')
  close all
  
% ENTER MODEL PARAMETERS
  repcost = 75;                                         % replacement cost  
  mancost =  5;                                         % maintanance cost
  delta = 0.9;                                          % discount factor  
  
% CONSTRUCT ACTION SPACE
  x = [0 1 2]';

% PACK MODEL STRUCTURE
  clear model  
  model.func = 'mfdp06';                                % model functions
  model.discount = delta;                               % discount factor
  model.actions = x;                                    % model actions
  model.discretestates = [1 2];                         % index of discrete states    
  model.params = {repcost mancost};                     % other parameters

% DEFINE APPROXIMATION SPACE
  fspace = fundefn([],[],[],[],[],[1:10]',[0:10]');
  scoord = funnode(fspace);	                            % state collocation grid coordinates
  s = gridmake(scoord);			                        % state collocation grid points
  
% SOLVE BELLMAN EQUATION
  [c,s,v,x,resid] = dpsolve(model,fspace,s);
  
% DISPLAY OPTIMAL REPLACEMENT-MAINTANENCE POLICY
  x1=reshape(x,10,11);
  fprintf('                       Optimal Replacement-Maintanence Policy\n')   
  fprintf('                                Number of Servicings\n')   
  fprintf('   Age %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i\n',scoord{2}');
  disp([scoord{1},x1])

