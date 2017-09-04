%% DEMDDP10 Stochastic Cow Replacement Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
delta = 0.9;                            % discount factor
cost  = 500;                            % replacement cost
price = 150;                            % milk price

% State space
s1 = (1:10)';                           % lactation states
s2 = [0.8;1.0;1.2];                     % productivity states
n1 = length(s1);                        % number of lactation states
n2 = length(s2);                        % number of productivity states
[S1,S2] = gridmake(s1,s2);              % combined state grid
n = n1*n2;                              % total number of states

% Action space (keep='K', replace='R')
X = ['K','R'];                          % keep or replace
m = length(X);                          % number of actions

% Reward function
y = (-0.2*S1.^2+2*S1+8).*S2;            % yield per lactation
f = [price*y price*y-cost];             % net revenue by action
f(S1==10,1) = -inf;                     % force replace at lactation 10

% State transition probability matrix
P = zeros(2,n1,n2,n1,n2);
for i=1:n1
  for j=1:n2
    if i<10
      P(1,i,j,i+1,j) = 1;     % Raise lactation number by 1, if keep
    else
      P(1,i,j,1,1) = 0.2;     % Forced replacement after lactation 10
      P(1,i,j,1,2) = 0.6;
      P(1,i,j,1,3) = 0.2;
    end
    P(2,i,j,1,1) = 0.2;       % Optional replacement
    P(2,i,j,1,2) = 0.6;
    P(2,i,j,1,3) = 0.2;
  end
end
P = reshape(P,2,n,n);

% Model structure
clear model
model.reward     = f;
model.transprob  = P;
model.discount   = delta;


%% SOLUTION

% Solve Bellman equation
[v,x,pstar] = ddpsolve(model);
   

%% ANALYSIS

% Display optimal policy
fprintf('\n\n')
fprintf('              Optimal Policy\n')
fprintf('     Age       Lo       Med        Hi\n')
fprintf('%8i %8c  %8c  %8c\n',[s1 reshape(X(x),n1,n2)]')

% Plot value function
figure
plot(s1,reshape(v,n1,n2))
legend('Low','Med','Hi')
title('Value Function')
xlabel('Age')
ylabel('Value of Dairy Operation')

% Steady-state distribution
P = markov(pstar);
fprintf('\n\n')
fprintf('          Invariant Distribution     \n')
fprintf('     Age       Lo       Med        Hi\n')
fprintf('%8i %8.3f  %8.3f  %8.3f\n',[s1 reshape(P,n1,n2)]')

% Steady-state mean cow age and productivity
avgage = P'*S1;
avgpri = P'*S2;
fprintf('\n\n')
fprintf('Steady-state Age          %8.2f\n',avgage)
fprintf('Steady-state Productivity %8.2f\n',avgpri)


%% SAVE FIGURES
printfigures(mfilename)