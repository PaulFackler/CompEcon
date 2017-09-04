%% DEMDOC05 Deterministic Production Adjustment Model
%
% Profit maximizing firm must decide how rapidly to adjust production.
%
% State
%     q       current production rate
% Control
%     x       production adjustment rate
% Parameters
%     alpha   production cost constant
%     beta    production cost elasticity
%     gamma   adjustment cost parameter
%     p       price
%     rho     continuous discount rate

function demdoc05

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha = 1;                                      % production cost constant
beta  = 1.5;                                    % production cost elasticity
gamma = 4;                                      % adjustment cost parameter
p     = 1;                                      % price
rho   = 0.1;                                    % continuous discount rate

% Model structure
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {p,alpha,beta,gamma};            % function file parameters

% Approximation structure
n = 10;                                         % number of collocation nodes
qmin = 0.2;                                     % minimum state
qmax = 0.7;                                     % maximum state
basis = fundefn('cheb',n,qmin,qmax);            % basis functions


%% SOLUTION

% Steady-state
qstar = (p/(beta*alpha))^(1/(beta-1));
xstar = 0;
vstar = (p*qstar-alpha*qstar.^beta)/rho;
fprintf('Steady States\n') 
fprintf('   Production    %5.2f\n'  ,qstar)

% Solve HJB equation by collocation
[c,q,v,x,resid] = docsolve(model,basis);


%% ANALYSIS

% Plot optimal policy
figure
hold on
plot(q,x)
title('Optimal Adjustment Policy')
xlabel('Production Rate')
ylabel('Adjustment Rate')

%  ... plot steady-state action
plothdash([],xstar)
plotvdash(qstar,xstar)
plotbullet(qstar,xstar)
plottext(qstar+0.01,[],'$q^*$')

% Plot value function
figure
hold on
plot(q,v)
title('Value Function')
xlabel('Production Rate')
ylabel('Value of the Firm')

%  ... plot steady-state value
plotvdash(qstar,vstar)
plotbullet(qstar,vstar)
plottext(qstar+0.01,[],'$q^*$')

% Plot shadow price function
figure
hold on
plot(q,funeval(c,basis,q,1))
title('Shadow Price Function')
xlabel('Production Rate')
ylabel('Shadow Price')

%  ... plot steady-state shadow price
pstar = funeval(c,basis,qstar,1);
plothdash(qstar,pstar)
plotvdash(qstar,pstar)
plotbullet(qstar,pstar)
plottext(qstar+0.01,[],'$q^*$')
plottext([],pstar,'$\lambda^*$')

% Plot residual
figure
hold on
plot(q,resid)
plothdash([],0)
title('HJB Equation Residual')
xlabel('Production Rate')
ylabel('Residual')


%% SIMULATION

% Initial state and time horizon
q0 = 0.2;            % initial production
T  = 10;             % time horizon

% Simulate model
[t,qsim] = docsimul(model,basis,q,x,q0,T);

% Time to midway adjustment
qhalf = (q0+qstar)/2;
thalf = interp1(qsim,t,qhalf,'pchip');

% Plot simulated state path
figure
hold on
plot(t,qsim)
title('Simulated Production')
xlabel('Time')
ylabel('Production')

%  ... plot steady-state
plothdash([],qstar)
plothdash(thalf,qhalf)
plotvdash(thalf,qhalf)
plotbullet(thalf,qhalf)
ylim([0.2 0.5])
plottext(0.08,qstar,'$q^*$')
plottext(0.08,qhalf,'$q_{0.5}$')
plottext(thalf+0.08,[],'$t_{0.5}$')


%% SAVE FIGURES
printfigures(mfilename)


%% DOCSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, and
%   transition at an arbitrary number of ns states:
%      out = func(flag,s,x,Vs,<params>)
%   where s is ns.ds states, x is ns.dx controls, Vs is ns.ds first
%   derivatives of value function.

function out = func(flag,q,x,Vq,p,alpha,beta,gamma)
switch flag
  case 'x'      % optimal control
    out = Vq/gamma;
  case 'f'      % reward
    k = alpha*q.^beta;
    a = 0.5*gamma*x.^2;
    out = p*q - k - a;
  case 'g'      % transition
    out = x;
end