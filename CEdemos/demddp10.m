% DEMDDP10 Stochastic Cow Replacement Model
  fprintf('\nDEMDDP10 STOCHASTIC COW REPLACEMENT MODEL\n')
  close all

% Enter model parameters
  delta = 0.9;                  % discount factor  
  cost  = 500;                  % replacement cost  
  price = 150;                  % milk price

% Construct state space
  s1 = (1:10)';                 % lactation states
  s2 = [0.8;1.0;1.2];           % productivity states
  n1 = length(s1);
  n2 = length(s2);
  [S1,S2] = gridmake(s1,s2);    % combined state grid
  n = n1*n2;                    % number of states

% Construct action space
  X = ['K','R'];                % keep or replace

% Construct reward function (actions keep=1, replace=2)
  y = (-0.2*S1.^2+2*S1+8).*S2;  % yield per lactation  
  f = [price*y price*y-cost];   % net revenue by action
  f(S1==10,1) = -inf;           % force replace at lactation 10

% Construct state transition probability matrix
  P = zeros(2,n1,n2,n1,n2);
  for i=1:n1
  for j=1:n2
     if i<10
        P(1,i,j,i+1,j) = 1;      % Raise lactation number by 1, if keep
     else
        P(1,i,j,1,1) = 0.2;      % Forced replacement after lactation 10
        P(1,i,j,1,2) = 0.6;
        P(1,i,j,1,3) = 0.2;
     end
     P(2,i,j,1,1) = 0.2;         % Optional replacement
     P(2,i,j,1,2) = 0.6;
     P(2,i,j,1,3) = 0.2;
  end
  end
  P = reshape(P,2,n,n);

% Pach model structure
  clear model
  model.reward     = f;
  model.transprob  = P;
  model.discount   = delta;

% Solve infinite-horizon model using policy iteration
  [v,x,pstar] = ddpsolve(model);

% Display optimal value
  disp('Value Function')
  disp('     Age       Lo       Med        Hi')
  fprintf('%8i %8.1f  %8.1f  %8.1f\n',[s1 reshape(v,n1,n2)]')

% Display optimal policy
  disp('Optimal Policy')
  disp('     Age       Lo       Med        Hi')
  fprintf('%8i %8c  %8c  %8c\n',[s1 reshape(X(x),n1,n2)]')

% Plot optimal value
  figure(1); plot(s1,reshape(v,n1,n2))
  xlabel('Age'); ylabel('Optimal Value');
  legend('Low','Med','Hi')

% Compute steady state distribution
  pi = markov(pstar);

% Display invariant distribution
  disp('          Invariant Distribution     ')
  disp('     Age       Lo       Med        Hi')
  fprintf('%8i %8.3f  %8.3f  %8.3f\n',[s1 reshape(pi,n1,n2)]')

% Compute steady-state cow age and productivity
  avgage = pi'*S1;
  avgpri = pi'*S2;
  fprintf('\nSteady-state Age          %8.2f\n',avgage)
  fprintf('\nSteady-state Productivity %8.2f\n',avgpri)
  
% Save Plots as EPS Files
  prtfigs(mfilename)
