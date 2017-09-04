%% DEMIC01 Asset Replacement Model

function demic01

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
beta = [1; 0.05; -0.003];
P    = 2;
C    = 3;
rho  = 0.1;

% Model structure
model.func   = @func;
model.params = {beta,P,rho};
model.xindex = [0 0;2 0];
model.F      = [0;C];


%% SOLUTION

% Initial values
x = [0 0;100 0];

% Solve collocation equation
n = 15;
[cv,basis,x] = icsolve(model,x,n);
Astar = x(2,1);


%% ANALYSIS

% Value function
A = nodeunif(101,0,Astar);
V = funeval(cv,basis,A);

% Plot results
figure
hold on
plot(A,V,'k')
title('Value Function')
xlabel('$A$')
ylabel('$V$')
xlim([0 18])


%% SAVE FIGURES
printfigures(mfilename)


%% ICSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, state
%   drift, and state diffusion at an arbitrary number of ns states:
%      out = func(flag,s,<params>)

function out = func(flag,s,beta,P,rho)
switch flag
  case 'f'        % reward
    Q = (beta(3)*s+beta(2)).*s+beta(1);
    out = P*Q;
  case 'mu'       % state drift
    out = ones(size(s,1),1);
  case 'sigma'    % state diffusion
    out = [];
  case 'rho'      % state contingent discount rate
    out = rho+zeros(size(s,1),1);
  case 'R+'       % reward associated with trigger s(1) and target s(2) (s(1)<s(2))
    out = zeros(1,4);
  case 'R-'       % reward associated with trigger s(1) and target s(2) (s(1)>s(2))
    out = zeros(1,4);
end