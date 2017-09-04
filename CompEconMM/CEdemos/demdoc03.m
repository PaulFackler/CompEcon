%% DEMDOC03 Deterministic Nonrenewable Resource Model
%
% Welfare maximizing social planner must decide the rate at which a
% nonrenewable resource should be harvested.
%
% State
%     s       resource stock
% Control
%     q       harvest rate
% Parameters
%     kappa   harvest unit cost scale factor
%     gamma   harvest unit cost elasticity
%     eta     inverse elasticity of demand
%     rho     continuous discount rate

function demdoc03

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
kappa = 10;                                     % harvest unit cost scale factor
gamma = 1;                                      % harvest unit cost elasticity
eta   = 1.5;                                    % inverse elasticity of demand
rho   = 0.05;                                   % continuous discount rate

% Model structure
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {eta,kappa,gamma};               % function file parameters

% Approximation structure
n = 20;                                         % number of basis functions
smin = 0.1;                                     % minimum state
smax = 1.0;                                     % maximum state
basis = fundefn('cheb',n,smin,smax);            % basis functions


%% SOLUTION

% Solve HJB equation by collocation
[c,s,v,q,resid] = docsolve(model,basis);


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
p = funeval(c,basis,s,1);
plot(s,p)
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

% Initial state and time horizon
s0 = smax;            % initial resource stock
T  = 40;              % time horizon

% Simulate Model
[t,ssim,qsim] = docsimul(model,basis,s,q,s0,T);

% Plot simulated state and control paths
figure
plot(t,[ssim qsim])
legend('Resource Stock','Rate of Harvest')
title('Simulated Resource Stock and Rate of Harvest')
xlabel('Time')
ylabel('Quantity')


%% SAVE FIGURES
printfigures(mfilename)


%% Model Function File
function out = func(flag,s,q,Vs,eta,kappa,gamma)
switch flag
  case 'x'      % optimal control
    k = kappa*s.^(-gamma);
    out = (Vs+k).^(-1/eta);
  case 'f'      % reward
    u = (1/(1-eta))*q.^(1-eta);
    k = kappa*s.^(-gamma);
    out = u - k.*q;
  case 'g'      % transition
    out = -q;
end