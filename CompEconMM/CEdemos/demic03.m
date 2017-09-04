%% DEMIC03 Storage Management Model

function demic03

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
mu    = -0.2;
sigma = 0.05;
k     = 0.05;
P     = 2;
F     = 3;
rho   = 0.1;

% Model structure
model.func   = @func;
model.params = {mu,sigma,k,P,rho};
model.xindex = [1 1;0 0];
model.F      = [F;0];


%% SOLUTION

% Initial values
smax = 8;
x = [0 smax;smax smax];

% Solve collocation equation
n = 10;
[cv,basis,x] = icsolve(model,x,n);
Sstar = x(1,2);


%% ANALYSIS

% Value function and derivatives
S = nodeunif(101,0,smax);
V = funeval(cv,basis,S);
Vstar = funeval(cv,basis,Sstar);

% Plot value function
figure
hold on
plot(S,V,S,Vstar+P*(S-Sstar))
legend('Chebychev Collocation','L-Q Approximation','Location','NW')
title('Value Function')
xlabel('S')
ylabel('V')

%  ... plot switch points
plotvdash(Sstar,Vstar)
plotbullet(Sstar,Vstar)
plottext(Sstar+0.01,[],'$S^*$')

% Display selected values
fprintf('Switching Points\n') 
fprintf('         S*    V(S*)\n')
fprintf('      %5.3f   %5.3f\n',Sstar,Vstar)


%% SAVE FIGURES
printfigures(mfilename)


%% ICSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, state
%   drift, and state diffusion at an arbitrary number of ns states:
%      out = func(flag,s,<params>)

function out = func(flag,s,mu,sigma,k,P,rho)
switch flag
  case 'f'        % reward
    out = -k*s;
  case 'mu'       % state drift
    out = mu+zeros(size(s,1),1);
  case 'sigma'    % state diffusion
    out = sigma+zeros(size(s,1),1);
  case 'rho'      % state contingent discount rate
    out = rho+zeros(size(s,1),1);
  case 'R+'       % reward associated with trigger s(1) and target s(2) (s(1)<s(2))
    out = [P*(s(1)-s(2)) P -P 0];
  case 'R-'       % reward associated with trigger s(1) and target s(2) (s(1)>s(2))
    out = zeros(1,4);
end