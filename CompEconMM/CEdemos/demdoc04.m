%% DEMDOC04 Deterministic Renewable Resource Model
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
%     kappa   unit harvest cost scale factor
%     gamma   unit harvest cost elasticity
%     eta     inverse elasticity of demand
%     rho     continuous discount rate

function demdoc04

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha = 0.25;                                   % biological growth function scale factor
beta  = 0.50;                                   % biological growth function elasticity
kappa = 5;                                      % unit harvest cost scale factor
gamma = 1.5;                                    % unit harvest cost elasticity
eta   = 1.5;                                    % inverse elasticity of demand
rho   = 0.05;                                   % continuous discount rate

% Model structure
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {alpha,beta,kappa,gamma,eta};    % function file parameters

% Approximation structure
n = 20;                                         % number of basis functions
smin = 0.2;                                     % minimum state
smax = 1.0;                                     % maximum state
basis = fundefn('cheb',n,smin,smax);            % basis functions


%% SOLUTION

% Steady-state
sstar = 0.6;
qstar = alpha*sstar.*(1-sstar.^beta);
xstar = [sstar;qstar];
xstar = broyden(@funcss,xstar,alpha,beta,kappa,gamma,eta,rho);
sstar = xstar(1);
qstar = xstar(2);
fprintf('Steady States\n') 
fprintf('   Resource Stock     %5.2f\n',sstar)
fprintf('   Rate of Harvest    %5.2f\n',qstar)

% Solve HJB equation by collocation
[c,s,v,q,resid] = docsolve(model,basis);


%% ANALYSIS

% Plot optimal policy
figure
hold on
plot(s,q)
title('Optimal Harvest Policy')
xlabel('Resource Stock')
ylabel('Rate of Harvest')

%  ... plot steady-state action
plothdash(sstar,qstar)
plotvdash(sstar,qstar)
plotbullet(sstar,qstar)
xl = xlim;
plottext(sstar+0.01,[],'$s^*$')
plottext(xl(1)+0.01,qstar,'$q^*$')

% Plot value function
figure
hold on
plot(s,v)
title('Value Function')
xlabel('Resource Stock')
ylabel('Social Welfare')

%  ... plot steady-state value
vst = funeval(c,basis,sstar);
plotvdash(sstar,vst)
plotbullet(sstar,vst)
plottext(sstar+0.01,[],'$s^*$')

% Plot shadow price function
figure
hold on
p = funeval(c,basis,s,1);
plot(s,p)
title('Shadow Price Function')
xlabel('Resource Stock')
ylabel('Shadow Price')

%  ... plot steady-state shadow price
pstar = funeval(c,basis,sstar,1);
plothdash(sstar,pstar)
plotvdash(sstar,pstar)
plotbullet(sstar,pstar)
plottext(sstar+0.01,[],'$s^*$')
plottext([],pstar,'$\lambda^*$')

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
s0 = smax;            % initial capital stock
T   = 20;             % time horizon

% Simulate model
[t,ssim,qsim] = docsimul(model,basis,s,q,s0,T);

% Plot simulated state and control paths
figure
hold on
plot(t,[ssim qsim])
legend('Resource Stock','Rate of Harvest')
title('Simulated Resource Stock and Rate of Harvest')
xlabel('Time')
ylabel('Quantity')

%  ... plot steady-state
plothdash([],sstar,'b')
plothdash([],qstar,'r')
plottext([],sstar,'$s^*$')
plottext([],qstar,'$q^*$')


%% SAVE FIGURES
printfigures(mfilename)


%% STEADY-STATE FILE

function out = funcss(x,alpha,beta,kappa,gamma,eta,rho)

s    = x(1,:);                             % map x to s
q    = x(2,:);                             % map x to q
p    = q.^(-eta);                          % inverse demand function
pder = -eta*q.^(-eta-1);                   % inverse demand derivative
g    = alpha*s.*(1-s.^beta);               % biological growth function
gder = alpha*(1-(1+beta)*s.^beta);         % marginal biological growth function
k    = kappa*s.^(-gamma);                  % harvest unit cost function
kder = -kappa*gamma*s.^(-gamma-1);         % harvest unit cost derivative
out  = [g-q; ((rho-gder).*(p-k)+kder.*g)./pder];


%% DOCSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, and
%   transition at an arbitrary number of ns states:
%      out = func(flag,s,x,Vs,<params>)
%   where s is ns.ds states, x is ns.dx controls, Vs is ns.ds first
%   derivatives of value function.

function out = func(flag,s,q,Vs,alpha,beta,kappa,gamma,eta)
switch flag
  case 'x'      % optimal control
    k = kappa*s.^(-gamma);
    out = (Vs+k).^(-1/eta);
  case 'f'      % reward
    u = (1/(1-eta))*q.^(1-eta);
    k = kappa*s.^(-gamma);
    out = u - k.*q;
  case 'g'      % transition
    g = alpha*s.*(1-s.^beta);
    out = g - q;
end