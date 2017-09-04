%% DEMDP09 Private Non-Renewable Resource Model
%
% Profit maximizing mine owner must decide how much ore to extract.
%
% States
%     s       ore stock
% Actions
%     q       ore extracted and sold
% Parameters
%     a       demand function parameters
%     b       cost function parameters
%     delta   discount factor

function demdp09

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
 
% Model parameters
a = [5 0.8];                            % demand function parameters
b = [7 1.0];                            % cost function parameters
delta = 0.9;                            % discount factor

% Model structure
model.func = @func;                     % model functions
model.params = {a b};                   % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions

% Approximation structure
n    = 101;                             % number of collocation nodes
smin =   0;                             % minimum state
smax =  10;                             % maximum state
basis = fundefn('spli',n,smin,smax);    % basis functions

% Check model derivatives
dpcheck(model,smax,0)


%% SOLUTION

% Solve collocation equation
[c,s,v,q,resid] = dpsolve(model,basis);

% Abandonment point
sstar = (b(1)-a(1))/b(2);
fprintf('Abandonment Point %5.2f\n',sstar)


%% ANALYSIS

% Plot optimal policy
figure
hold on
plot(s,q)
title('Optimal Extraction Policy')
xlabel('Ore Stock')
ylabel('Quantity Extracted')

%  ... plot abandonment level
plotvdash(sstar,0)
plotbullet(sstar,0)
plottext(sstar,[],'$s^*$')

% Plot value function
figure
plot(s,v)
title('Value Function')
xlabel('Ore Stock')
ylabel('Value of Mine')

% Plot shadow price function
figure
hold on
p = funeval(c,basis,s,1);
plot(s,p)
title('Shadow Price Function')
xlabel('Ore Stock')
ylabel('Shadow Price')

%  ... plot abandonment level
plotvdash(sstar,0)
plotbullet(sstar,0)
plottext(sstar,[],'$s^*$')

% Plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
title('Bellman Equation Residual')
xlabel('Ore Stock')
ylabel('Residual')

%  ... plot abandonment level
plotvdash(sstar,0)
plotbullet(sstar,0)
plottext(sstar,[],'$s^*$')


%% SIMULATION

% Simulation parameters
nper = 21;                              % number of periods simulated
nrep = 	1;                              % number of replications

% Initialize simulation
sinit = smax;

% Simulate model
[ssim,qsim] = dpsimul(model,basis,nper,sinit,[],s,v,q);

% Plot simulated state and policy path
figure
hold on
plot(0:nper-1,ssim,0:nper-1,qsim)
legend('Ore Stock','Quantity Extracted')
title('Simulated Ore Stock and Quantity Extracted')
xlabel('Period')
ylabel('Quantity')

%  ... plot abandonment level
plothdash([],sstar)
plottext([],sstar+0.01,'$s^*$')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions according to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,e,<params>)
%    where s is ns.ds continuous states, x is ns.dx continuous actions, i
%    is ns.1 or scalar discrete states, j is ns.1 or scalar discrete
%    actions, and e is ns.de continuous state transition shocks.

function [out1,out2,out3] = func(flag,s,q,i,j,e,a,b)

switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = s;
    out3 = [];
  case 'f'      % reward
    out1 = (a(1)-b(1)+b(2)*s).*q - (a(2)+b(2)/2).*q.^2;
    out2 = (a(1)-b(1)+b(2)*s) - 2*(a(2)+b(2)/2).*q;
    out3 = -2*(a(2)+b(2)/2)*ones(size(s));
  case 'g'      % transition
    out1 = s-q;
    out2 = -ones(size(s));
    out3 = zeros(size(s));
end