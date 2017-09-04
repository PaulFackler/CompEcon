%% DEMIC02 Timber Harvesting Model

function demic02

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
alpha = 0.1;
m     = 1;
sigma = 0.05;
P     = 3;
C     = 0.15;
rho   = 0.1;

% Model structure
model.func   = @func;
model.params = {alpha,m,sigma,P,rho};
model.xindex = [0 0;2 0];
model.F      = [0;C];


%% SOLUTION

% Initial values
x = [0 0;0.5 0];

% Solve collocation equation
n = 15;
[cv,basis,x] = icsolve(model,x,n);
Sstar = x(2,1);


%% ANALYSIS

% Value function and derivatives
S = nodeunif(101,0,m/2);
V = funeval(cv,basis,S);
Vstar = funeval(cv,basis,Sstar);
dV = funeval(cv,basis,S,1);
dVstar = funeval(cv,basis,Sstar,1);
V(S>Sstar) = Vstar+P*(S(S>Sstar)-Sstar);
dV(S>Sstar) = P;

% Plot value function
figure
hold on
plot(S,V)
title('Value Function')
xlabel('$S$')
ylabel('$V$')
ylim([1.8 3.2])

%  ... plot switch points
plotvdash(Sstar,Vstar)
plotbullet(Sstar,Vstar)
plottext(Sstar+0.01,[],'$S^*$')

% Plot first derivative of value function
figure
hold on
plot(S,dV)
title('First Derivative of Value Function')
xlabel('$S$')
ylabel('$V''$')
ylim([1.8 3.2])

%  ... plot switch points
plotvdash(Sstar,dVstar)
plotbullet(Sstar,dVstar)
plottext(Sstar+0.01,[],'$S^*$')


%% SAVE FIGURES
printfigures(mfilename)


%% ICSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, state
%   drift, and state diffusion at an arbitrary number of ns states:
%      out = func(flag,s,<params>)

function out = func(flag,s,alpha,m,sigma,P,rho)
switch flag
  case 'f'        % reward
    out = zeros(size(s,1),1);
  case 'mu'       % state drift
    out = alpha*(m-s);
  case 'sigma'    % state diffusion
    out = sigma*sqrt(s);
  case 'rho'      % state contingent discount rate
    out = rho+zeros(size(s,1),1);
  case 'R+'       % reward associated with trigger s(1) and target s(2) (s(1)<s(2))
    out = zeros(1,4);
  case 'R-'       % reward associated with trigger s(1) and target s(2) (s(1)>s(2))
    out = [P*(s(1)-s(2)) P -P 0];
end