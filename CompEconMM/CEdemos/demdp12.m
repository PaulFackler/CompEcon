%% DEMDP12 Production Management Model
%
% Profit maximizing entrepeneur must decide how much to produce, subject to 
% production adjustment costs.
%
% States
%     i       market price (discrete)
%     s       lagged production (continuous)
% Actions
%     x       current production
% Parameters
%     alpha   marginal adjustment cost
%     beta    marginal production cost parameters
%     pbar    long-run average market price
%     sigma   market price shock standard deviation
%     delta   discount factor

function demdp12

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha = 0.01;                               % marginal adjustment cost
beta  = [0.8 0.03];                         % marginal production cost parameters
pbar  = 1.0;                                % long-run average market price
sigma = 0.2;                                % market price shock standard deviation
delta = 0.9;                                % discount factor

% Continuous state shock distribution
m = 3;                                      % number of shocks
mu = log(pbar)-sigma^2/2;                   % mean log price
[p,w] = qnwlogn(m,mu,sigma^2);              % shocks and probabilities
q = w(:,ones(1,m))';                        % transition probabilities

% Model structure
model.func = @func;                         % model functions
model.params = {alpha beta p};              % function parameters
model.discount = delta;                     % discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 1;                               % dimension of continuous action
model.ni = m;                               % number of discrete states
model.nj = 0;                               % number of discrete actions
model.q  = q;                               % discrete state transition probabilities

% Approximation structure
n = 50;                                     % number of collocation nodes
smin =  0;                                  % minimum state
smax = 20;                                  % maximum state
basis = fundefn('spli',n,smin,smax);        % basis functions


%% SOLUTION

% Deterministic steady-state
sstar = (pbar-beta(1))/beta(2);

% Check model derivatives
dpcheck(model,sstar,sstar)

% Solve collocation equation
[c,s,v,x,resid] = dpsolve(model,basis);


%% ANALYSIS

% Plot optimal policy
figure
plot(s,x)
legend('Low Price','Average Price','High Price')
title('Optimal Production Policy')
xlabel('Lagged Production')
ylabel('Production')

% Plot value function
figure
plot(s,v)
legend('Low Price','Average Price','High Price')
title('Value Function')
xlabel('Lagged Production')
ylabel('Value of the Firm')

% Plot shadow price function
figure
lambda = funeval(c,basis,s,1);
plot(s,lambda)
legend('Low Price','Average Price','High Price')
title('Shadow Price of Lagged Production')
xlabel('Lagged Production')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
legend('Low Price','Average Price','High Price')
title('Bellman Equation Residual')
xlabel('Lagged Production')
ylabel('Residual')


%% SIMULATION

% Simulation parameters
nper = 26;                              % number of periods simulated
nrep = 10000;                           % number of replications

% Initialize simulation
sinit = smin(ones(nrep,1),:);           % initial lagged production
iinit = 2*ones(nrep,1);                 % initial market price state
rng('default')

% Simulate model
[ssim,xsim,isim] = dpsimul(model,basis,nper,sinit,iinit,s,v,x);
psim = p(isim);

% Ergodic moments
savg = mean(ssim(:)); 
pavg = mean(psim(:)); 
sstd = std(ssim(:)); 
pstd = std(psim(:)); 
fprintf('Ergodic Moments\n') 
fprintf('          Deterministic    Ergodic      Ergodic\n') 
fprintf('          Steady-State      Mean     Std Deviation\n') 
fprintf('Price        %5.3f         %5.3f         %5.3f\n'  ,[pbar  pavg pstd])
fprintf('Production   %5.3f         %5.3f         %5.3f\n\n',[sstar savg sstd])

% Plot simulated action path
figure
hold on
plot(0:nper-1,xsim(1:3,:))
plot(0:nper-1,mean(xsim),'k')
title('Simulated and Expected Production')
xlabel('Period')
ylabel('Production')


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

function [out1,out2,out3] = func(flag,s,q,i,j,e,alpha,beta,p)

n = length(s);
p = p(i);
l = s;
switch flag
  case 'b'      % bounds
    out1 = zeros(n,1);
    out2 = inf*ones(n,1);
    out3 = [];
  case 'f'      % reward
    out1 = p.*q - (beta(1)*q+0.5*beta(2)*q.^2) - 0.5*alpha*((q-l).^2);
    out2 = p - beta(1) - beta(2)*q - alpha*(q-l);
    out3 = (-beta(2)-alpha)*ones(n,1);
  case 'g'      % transition
    out1 = q;
    out2 = ones(n,1);
    out3 = zeros(n,1);
end