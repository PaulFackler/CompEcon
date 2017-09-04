%% DEMSOC06 Stochastic Production-Adjustment Model
%
% Profit maximizing entrepreneur must decide how to adjust production,
% given stochastic prices and adjustment costs.
%
% State
%     p       output price
%     q       production rate
% Control
%     x       production adjustment rate
% Parameters
%     alpha   production cost constant
%     beta    production cost elasticity
%     gamma   adjustment cost parameter
%     pbar    mean price
%     delta   speed of price mean reversion
%     sigma   price volatility
%     rho     discount rate

function demsoc06

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha = 1;                                      % production cost constant
beta  = 1.5;                                    % production cost elasticity
gamma = 4;                                      % adjustment cost parameter
pbar  = 1;                                      % mean price
delta = 0.5;                                    % speed of price mean reversion
sigma = 0.2;                                    % price volatility
rho   = 0.1;                                    % discount rate

% Model structure
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {delta,pbar,sigma,gamma,alpha,beta};  % function file parameters

% Approximation structure
n = [15 11];                                    % number of basis functions per state dimension
pmin = 0.1;                                     % minimum output price        
pmax = 1.6;                                     % maximum output price
qmin = 0.1;                                     % minimum production rate
qmax = 1.0;                                     % maximum production rate
smin = [pmin qmin];                             % minimum states
smax = [pmax qmax];                             % maximum states
basis = fundefn('cheb',n,smin,smax);            % basis functions


%% SOLUTION

% Solve HJB equation by collocation
[c,s,v,x,resid] = socsolve(model,basis);

% Reshape output for plotting
nr = n*10+1;
p = reshape(s(:,1),nr);
q = reshape(s(:,2),nr);
v = reshape(v,nr);
x = reshape(x,nr);
resid = reshape(resid,nr);

% Shadow prices
p1 = funeval(c,basis,s,[1 0]);
p1 = reshape(p1,nr);
p2 = funeval(c,basis,s,[0 1]);
p2 = reshape(p2,nr);


%% ANALYSIS

% Plot optimal policy
figure
mesh(p,q,x)
title('Optimal Production Adjustment Policy')
xlabel('Output Price')
ylabel('Production Rate')
zlabel('Rate of Adjustment')

% Plot value function
figure
mesh(p,q,v)
title('Value Function')
xlabel('Output Price')
ylabel('Production Rate')
zlabel('Value of the Firm')

% Plot shadow price function 1
figure
mesh(p,q,p1)
title('Shadow Price of Output Price')
xlabel('Output Price')
ylabel('Production Rate')
zlabel('Shadow Price')

% Plot shadow price function 2
figure
mesh(p,q,p2)
title('Shadow Price of Production Level')
xlabel('Output Price')
ylabel('Production Rate')
zlabel('Shadow Price')

% Plot residual
figure
mesh(p,q,resid)
title('HJB Equation Residual')
xlabel('Output Price')
ylabel('Production Rate')
zlabel('Residual')


%% SIMULATION

% Initial state, time horizon, and number of replications
sinit = [pbar qmin];    % initial capital stock
T     = 25;             % time horizon
nrep = 5000;            % number of replications

% Simulate model
[t,ssim,xsim,smean,xmean] = socsimul(c,model,basis,sinit,T,nrep);
psim = ssim(:,:,1);
qsim = ssim(:,:,2);
pmean = smean(:,1);
qmean = smean(:,2);

% Plot simulated and expected state path
figure
plot(t,psim,t,pmean,'k')
title('Simulated and Expected Output Price')
xlabel('Time')
ylabel('Output Price')

% Plot simulated and expected state path
figure
plot(t,qsim,t,qmean,'k')
title('Simulated and Expected Production Rate')
xlabel('Time')
ylabel('Production Rate')

% Plot simulated and expected control path
figure
plot(t,xsim,t,xmean,'k')
title('Simulated and Expected Rate of Adjustment')
xlabel('Time')
ylabel('Rate of Adjustment')


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

function out = func(flag,s,x,Vs,~,delta,pbar,sigma,gamma,alpha,beta)
n = size(s,1);
p = s(:,1);
q = s(:,2);
switch flag
  case 'x'        % optimal control
    out = Vs(:,2)/gamma;
  case 'f'        % reward
    k = alpha*q.^beta;
    a = 0.5*gamma*x.^2;
    out = p.*q - k - a;
  case 'mu'       % state drift
    out = [delta*(pbar-p) x];
  case 'sigma'    % state diffusion
    out = zeros(n,2,2);
    out(:,1,1) = sigma*p;
end