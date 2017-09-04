%% DEMSOC04 Stochastic Renewable Resource Model
%
% Welfare maximizing social planner must decide the rate at which a
% renewable resource should be harvested.
%
% State
%     s       resource stock
% Control
%     q       harvest rate
% Parameters
%     alpha   biological growth function scale factor
%     beta    biological growth function elasticity
%     kappa   harvest unit cost scale factor
%     gamma   harvest unit cost elasticity
%     eta     inverse elasticity of demand
%     sigma   biological growth volatility
%     rho     discount rate

function demsoc04

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha = 0.5;                                    % biological growth function scale factor
beta  = 0.6;                                    % biological growth function elasticity
kappa = 5.0;                                    % harvest unit cost scale factor
gamma = 2.0;                                    % harvest unit cost elasticity
eta   = 2.0;                                    % inverse elasticity of demand
sigma = 0.05;                                   % biological growth volatility
rho   = 0.05;                                   % discount rate

% Model structure
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {alpha,beta,kappa,gamma,eta,sigma}; % function file parameters

% Approximation structure
n = 30;                                         % number of basis functions
smin = 0.2;                                     % minimum resource stock
smax = 1.0;                                     % maximum resource stock
basis = fundefn('cheb',n,smin,smax);            % basis functions


%% SOLUTION

% Solve HJB equation by collocation
[c,s,v,q,resid] = socsolve(model,basis);


%% ANALYSIS

% Plot optimal policy
figure
plot(s,q)
title('Optimal Harvest Policy')
xlabel('Resource Stock')
ylabel('Rate of Harvest')

% Plot value function
figure
plot(s,v)
title('Value Function')
xlabel('Resource Stock')
ylabel('Social Welfare')

% Plot shadow price function
figure
Vs = funeval(c,basis,s,1);
plot(s,Vs)
title('Shadow Price Function')
xlabel('Resource Stock')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
title('HJB Equation Residual')
xlabel('Resource Stock')
ylabel('Residual')


%% SIMULATION

% Initial state, time horizon, and number of replications
sinit = smin;           % initial capital stock
T     = 60;             % time horizon
nrep  = 1000;           % number of replications

% Simulate model
[t,ssim,xsim,smean,xmean] = socsimul(c,model,basis,sinit,T,nrep);

% Plot simulated and expected state path
figure
plot(t,ssim,t,smean,'k')
title('Simulated and Expected Resource Stock')
xlabel('Time')
ylabel('Resource Stock')

% Plot simulated and expected control path
figure
plot(t,xsim,t,xmean,'k')
title('Simulated and Expected Rate of Harvest')
xlabel('Time')
ylabel('Rate of Harvest')

% Plot ergodic distribution
figure
cp = itodensity(model,basis,c);
p  = max(0,funeval(cp,basis,s));
plot(s,p)
title('Ergodic Distribution of Resource Stock')
xlabel('Resource Stock')
ylabel('Probability')


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

function out = func(flag,s,q,Vs,~,alpha,beta,kappa,gamma,eta,sigma)
switch flag
  case 'x'        % optimal control
    k = kappa*s.^(-gamma);
    out = (Vs+k).^(-1/eta);
  case 'f'        % reward
    u = (1/(1-eta))*q.^(1-eta);
    k = kappa*s.^(-gamma);
    out = u - k.*q;
  case 'mu'       % state drift
    g = alpha*s.*(1-s.^beta);
    out = g - q;
  case 'sigma'    % state diffusion
    out = sigma*s;
end


%% PLOT ERRORS FOR CASES WITH ANALYTIC SOLUTION
% 
% % Plot Percent Approximation Error for Cases with Closed-Form Solution
% if gamma==(1+beta)&&eta==(1+beta)
%   optset('broyden','tol',1e-14)
%   optset('broyden','showiters',1)
%   theta = ((alpha+rho)*beta/(1+beta)-0.5*beta*sigma^2)^((1+beta)/beta);
%   fphi = @(phi) theta*phi^((1+beta)/beta)-beta*phi-kappa;
%   optset('broyden','tol',1e-12)
%   phi = broyden(fphi,12);
%   va = -phi*(s.^-beta+alpha/rho);
%   qa = ((kappa+beta*phi)^(-1/(1+beta)))*s;
%   figure
%   hold on
%   plot(s,v./va-1,s,q./qa-1)
%   legend('Value Function','Optimal Policy')
%   plot(s,0*s,'k--','Linewidth',1)
%   title('Percent Approximation Errors')
%   xlabel('Resource Stock')
%   ylabel('Error')
% end