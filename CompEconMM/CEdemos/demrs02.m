%% DEMRS02 Optimal Fish Harvest Model
%
% Profit maximizing fisheries owner must decide harvesting effort to exert.
%
% State
%     s       stock of fish
% Control
%     h       harvest effort
% Parameters
%     alpha   biological growth function scale factor
%     sigma   biological growth volatility
%     H       maximum harvest effort
%     p       market price of fish
%     k       cost function parameter
%     rho     continuous discount rate

function demrs02

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
alpha = 0.5;                                    % biological growth function scale factor
sigma = 0.5;                                    % biological growth volatility
H     = 1;                                      % maximum harvest effort
p     = 1;                                      % market price of fish
k     = 0.25;                                   % cost function parameter
rho   = 0.1;                                    % continuous discount rate

% Model structure
model.func   = @func;
model.params = {alpha,sigma,H,p,k,rho};
model.xindex = [
  1 0 0 0 0;
  1 2 1 1 2;
  2 0 0 0 0];


%% SOLUTION

% Initial values and approximation size
x = [0.001;0.5;3];
n = [75 20];

% Solve collocation equation
[cv,basis,x] = rssolve(model,x,n);

sstar  = x(2);
vstar  = funeval(cv{1},basis{1},sstar);
vstar1 = funeval(cv{1},basis{1},sstar,1);
vstar2 = funeval(cv{1},basis{1},sstar,2);


%% ANALYSIS

% Value function and derivatives
s = nodeunif(2001,0.001,1);
m = 2;
Vi = zeros(size(s,1),m,3);
for i = 1:2
  Vi(:,i,1) = funeval(cv{i},basis{i},s);
  Vi(:,i,2) = funeval(cv{i},basis{i},s,1);
  Vi(:,i,3) = funeval(cv{i},basis{i},s,2);
end
V   = [Vi(s<sstar,1,1);Vi(s>=sstar,2,1)];
Vs  = [Vi(s<sstar,1,2);Vi(s>=sstar,2,2)];
Vss = [Vi(s<sstar,1,3);Vi(s>=sstar,2,3)];

% Plot value function
figure
hold on
plot(s,V)
title('Value Function')
xlabel('Fish Stock')
ylabel('$V$')

%  ... plot switch point
plotvdash(sstar,vstar)
plotbullet(sstar,vstar)
plottext(sstar+0.01,[],'$s^*$')

% Plot first derivative of value function
figure
hold on
plot(s,Vs)
ylim([0 5])
title('First Derivative of Value Function')
xlabel('Fish Stock')
ylabel('$V''$')

%  ... plot switch point
plotvdash(sstar,vstar1)
plotbullet(sstar,vstar1)
plottext(sstar+0.01,[],'$s^*$')

% Plot second derivative of value function
figure
hold on
plot(s,Vss)
title('First Derivative of Value Function')
xlabel('Fish Stock')
ylabel('$V''''$')
ylim([-10 5])

%  ... plot switch point
plotvdash(sstar,vstar2)
plotbullet(sstar,vstar2)
plottext(sstar+0.01,[],'$s^*$')

% Plot residual
e = rho*V-H*(p-k./s).*s.*(s>sstar)-(alpha*(1-s)-(s>sstar)*H).*s.*Vs-0.5*sigma^2*s.^2.*Vss;
figure
plot(s,e)
title('Approximation Residual')
xlabel('Fish Stock')
ylabel('Residual')

% print switching point
fprintf('\nSwitching point\n')
fprintf('  s*                 %5.3f\n',sstar)


%% SAVE FIGURES
printfigures(mfilename)


%% RSSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, state
%   drift, and state diffusion at an arbitrary number of ns states:
%      out = func(flag,s,x,<params>)

function out = func(flag,s,x,alpha,sigma,H,p,k,rho)

switch flag
  case 'f'        % reward
    out = (p-k./s).*s.*(x==2)*H;
  case 'g'        % state drift
    out = (alpha*(1-s)-(x==2)*H).*s;
  case 'sigma'    % state diffusion
    out = sigma.*s; 
  case 'rho'      % state contingent discount rate
    out = rho+zeros(size(s,1),1);
  case 'reward'   % jump rewards and derivatives
    out = [0 0 0;0 0 0;0 0 0];
end