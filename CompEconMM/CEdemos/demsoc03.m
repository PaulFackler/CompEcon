%% DEMSOC03 Stochastic Optimal Economic Growth Model
%
% Social benefit maximizing social planner must decide how much society
% should consume and invest.
%
% State
%     k       capital stock
%     y       productivity shock
% Control
%     q       consumption rate
% Parameters
%     alpha   capital share
%     delta   capital depreciation rate
%     theta   relative risk aversion
%     gamma   productivity mean reversion coefficient
%     sigma   productivity volatility
%     rho     discount rate

function demsoc03

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha =  0.4;                                   % capital share
delta =  0.1;                                   % capital depreciation rate
theta =  2.0;                                   % relative risk aversion
gamma =  0.5;                                   % productivity mean reversion coefficient
sigma =  0.05;                                  % productivity volatility
rho   =  0.05;                                  % discount rate

% Model Structure
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {alpha,delta,gamma,theta,sigma}; % function file parameters

% Approximation Structure
n = [15 5];                                     % number of basis functions per state dimension
kmin =  2;                                      % minimum capital stock
kmax =  8;                                      % maximum capital stock
ymin = 0.9;                                     % minimum productivity factor
ymax = 1.1;                                     % maximum productivity factor
smin = [kmin ymin];                             % minimum states
smax = [kmax ymax];                             % maximum states
[basis,~,snode] = fundefn('cheb',n,smin,smax);  % basis functions


%% SOLUTION

% Solve HJB equation by collocation
v = (((rho*snode(:,1)).^(1-theta))/(1-theta));
[c,s,v,q,res] = socsolve(model,basis,v);


%% ANALYSIS

% Fix y=1
p = funeval(c,basis,s,[1 0]);
k = s(:,1);
y = s(:,2);
j = find(y==1);
k = k(j);
v = v(j);
q = q(j);
p = p(j);
res = res(j);

% Plot optimal policy
figure
plot(k,q)
title('Optimal Consumption Policy')
xlabel('Capital Stock')
ylabel('Rate of Consumption')

% Plot value function
figure
plot(k,v)
title('Value Function')
xlabel('Capital Stock')
ylabel('Lifetime Utility')

% Plot shadow price function
figure
plot(k,p)
title('Shadow Price Function')
xlabel('Capital Stock')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(k,res)
plothdash([],0)
title('HJB Equation Residual')
xlabel('Capital Stock')
ylabel('Residual')


%% SIMULATION

% Initial state, time horizon, and number of replications
sinit = [3 1];          % initial capital stock and productivity shock
T     = 25;             % time horizon
nrep  = 2000;           % number of replications

% Simulate model
[t,ssim,xsim,smean,xmean] = socsimul(c,model,basis,sinit,T,nrep);
ksim = ssim(:,:,1);
ysim = ssim(:,:,2);
kmean = smean(:,1);
ymean = smean(:,2);

% Plot simulated and expected state path
figure
plot(t,ksim,t,kmean,'k')
title('Simulated and Expected Capital Stock')
xlabel('Time')
ylabel('Capital Stock')

% Plot simulated and expected state path
figure
plot(t,ysim,t,ymean,'k')
title('Simulated and Expected Productivity Shock')
xlabel('Time')
ylabel('Productivity Shock')

% Plot expected control path
figure
plot(t,xsim,t,xmean,'k')
title('Simulated and Expected Rate of Consumption')
xlabel('Time')
ylabel('Rate of Consumption')


%% SOLVE HJB EQUATION DIRECTLY USING BROYDEN

% Approximation Structure
n = [15 5];
smin = [ 6 0.9];
smax = [14 1.1];
[basis,Phi,s] = fundefn('cheb',n,smin,smax);

% Define Residual Function
k     = @(s) s(:,1);
y     = @(s) s(:,2);
V     = @(c,s) funeval(c,basis,s);
Vk    = @(c,s) funeval(c,basis,s,[1 0]);
Vy    = @(c,s) funeval(c,basis,s,[0 1]);
Vyy   = @(c,s) funeval(c,basis,s,[2 2]);
q     = @(c,s) Vk(c,s).^(-1/theta);
resid = @(c,s) (q(c,s).^(1-theta))/(1-theta) ...
       + (y(s).*k(s).^alpha-delta*k(s)-q(c,s)).*Vk(c,s) ...
       + gamma*(1-y(s)).*Vy(c,s) ...
       - 0.5*(sigma^2)*y(s).*Vyy(c,s) - rho*V(c,s);

% Solve
v = (((rho*k(s)).^(1-theta))/(1-theta));
c = Phi\v;
c = broyden(resid,c,s);
norm(resid(c,s))


%% SAVE FIGURES
printfigures(mfilename)


%% SOCSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, state
%   drift, and state diffusion at an arbitrary number of ns states:
%      out = func(flag,s,x,Vs,Vss,<params>)
%   where s is ns.ds states, x is ns.dx controls, Vs is ns.ds first
%   derivatives of value function, and Vss is ns.ds.ds second derivatives
%   of value function.

function out = func(flag,s,q,Vs,~,alpha,delta,gamma,theta,sigma)
k = s(:,1);
y = s(:,2);
switch flag
  case 'x'        % optimal control
    Vk  = Vs(:,1);
    out = Vk.^(-1/theta);
  case 'f'        % reward
    out = (1/(1-theta))*q.^(1-theta);
  case 'mu'       % state drift
    f = k.^alpha;
    out = [(y.*f-delta*k-q)  gamma*(1-y)];
  case 'sigma'    % state diffusion
    n = length(k);
    out = zeros(n,2,2);
    out(:,2,2) = sigma*sqrt(y);
end