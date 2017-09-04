%% DEMDP07 Stochastic Optimal Economic Growth Model
%
% Welfare maximizing social planner must decide how much society should
% consume and invest. Unlike the Ramsey model, this model allows arbitrary
% constant relative risk aversion, capital depreciation, and stochastic
% production shocks.  Also unlike the Ramsey model, it lacks a known
% closed-form solution.
%
% States
%     s       stock of wealth
% Actions
%     k       capital investment
% Parameters
%     alpha   relative risk aversion
%     beta    capital production elasticity
%     gamma   capital survival rate
%     sigma   production shock volatility
%     delta   discount factor

function demdp07

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
alpha = 0.2;                            % relative risk aversion
beta  = 0.5;                            % capital production elasticity
gamma = 0.9;                          	% capital survival rate
sigma = 0.1;                            % production shock volatility
delta = 0.9;                            % discount factor

% Continuous state shock distribution
m = 5;                                 	% number of shocks
[e,w] = qnwlogn(m,-sigma^2/2,sigma^2); 	% shocks and probabilities

% Model structure
model.func = @func;                     % model functions
model.params = {alpha beta gamma};	    % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities

% Approximation structure
n     = 10;                             % number of collocation nodes
smin  =  4;                             % minimum state
smax  = 12;                             % maximum state
basis = fundefn('cheb',n,smin,smax);    % basis functions


%% SOLUTION
  
% Deterministic steady-state
kstar = ((1-delta*gamma)/(delta*beta))^(1/(beta-1));  % capital investment
sstar = gamma*kstar + kstar^beta;                     % wealth

% Check model derivatives
dpcheck(model,sstar,kstar)

% Solve collocation equation
[c,s,v,k,resid] = dpsolve(model,basis);


%% ANALYSIS
 
% Plot optimal policy
figure
hold on
plot(s,k)
title('Optimal Investment Policy')
xlabel('Wealth')
ylabel('Capital Investment')

% Plot value function
figure
plot(s,v)
title('Value Function')
xlabel('Wealth')
ylabel('Lifetime Utility')

% Plot shadow price function
figure
pr = funeval(c,basis,s,1);
plot(s,pr)
title('Shadow Price Function')
xlabel('Wealth')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
title('Bellman Equation Residual')
xlabel('Wealth')
ylabel('Residual')


%% SIMULATION

% Simulation parameters
nper = 21;                              % number of periods simulated
nrep = 50000;                           % number of replications

% Initialize simulation
sinit = smin*ones(nrep,1);              % initial wealths
rng('default')

% Simulate model
[ssim,ksim] = dpsimul(model,basis,nper,sinit,[],s,v,k);

% Plot simulated and expected state path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Wealth')
xlabel('Period')
ylabel('Wealth')

% Plot simulated and expected action path
figure
hold on
plot(0:nper-1,ksim(1:3,:))
plot(0:nper-1,mean(ksim),'k')
title('Simulated and Expected Investment')
xlabel('Period')
ylabel('Investment')

% Ergodic moments
ssim = ssim(:,nper);
ksim = ksim(:,nper);
savg = mean(ssim); 
kavg = mean(ksim); 
sstd = std(ssim); 
kstd = std(ksim); 
fprintf('Ergodic Moments\n') 
fprintf('          Deterministic    Ergodic      Ergodic\n') 
fprintf('          Steady-State      Mean     Std Deviation\n') 
fprintf('Wealth       %5.3f         %5.3f         %5.3f\n'  ,[sstar savg sstd])
fprintf('Investment   %5.3f         %5.3f         %5.3f\n\n',[kstar kavg kstd])

% Plot ergodic state distribution
[qq,ss] = ksdensity(ssim(:),'support','positive');
figure
plot(ss,qq)
title('Ergodic Wealth Distribution')
xlabel('Wealth')
ylabel('Probability')
xlim([5 10])  


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

function [out1,out2,out3] = func(flag,s,k,~,~,e,alpha,beta,gamma)

switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = 0.99*s;
    out3 = [];
  case 'f'      % reward
    out1 = ((s-k).^(1-alpha))/(1-alpha);
    out2 = -(s-k).^(-alpha);
    out3 = -alpha*(s-k).^(-alpha-1);
  case 'g'      % transition
    out1 = gamma*k + e.*k.^beta;
    out2 = gamma + beta*e.*k.^(beta-1);
    out3 = (beta-1)*beta*e.*k.^(beta-2);
end