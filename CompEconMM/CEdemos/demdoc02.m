%% DEMDOC02 Deterministic Optimal Economic Growth Model
%
% Social benefit maximizing social planner must decide how much society
% should consume and invest.
%
% State
%     k       capital stock
% Control
%     q       consumption rate
% Parameters
%     alpha   capital share
%     delta   capital depreciation rate
%     theta   relative risk aversion
%     rho     continuous discount rate

function demdoc02

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha =  0.4;                                   % capital share
delta =  0.1;                                   % capital depreciation rate
theta =  2.0;                                   % relative risk aversion
rho   =  0.05;                                  % continuous discount rate

% Model structure
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {alpha,delta,theta};             % function file parameters

% Approximation structure
n = 21;                                         % number of basis functions 
kmin = 1;                                       % minimum state
kmax = 7;                                       % maximum state
[basis,~,k] = fundefn('cheb',n,kmin,kmax);      % basis functions


%% SOLUTION

% Steady-state
kstar = ((delta+rho)/alpha)^(1/(alpha-1));      % capital stock
qstar = kstar.^alpha-delta*kstar;               % consumption rate
vstar = ((1/(1-theta))*qstar.^(1-theta))/rho;   % value
lstar =  qstar.^(-theta);                       % shadow price
fprintf('Steady States\n') 
fprintf('   Capital Stock          %5.2f\n',kstar)
fprintf('   Rate of Consumption    %5.2f\n',qstar)
fprintf('   Shadow Price           %5.2f\n',lstar)

% Solve HJB equation by collocation
v = (((rho*k).^(1-theta))/(1-theta));
[c,k,v,q,res] = docsolve(model,basis,v);


%% ANALYSIS

% Plot optimal policy
figure
hold on
plot(k,q)
title('Optimal Consumption Policy')
xlabel('Capital Stock')
ylabel('Rate of Consumption')

%  ... plot steady-state action
plothdash(kstar,qstar)
plotvdash(kstar,qstar)
plotbullet(kstar,qstar)
plottext(kstar,[],'$k^*$')
plottext([],qstar,'$q^*$')

% Plot value function
figure
hold on
plot(k,v)
title('Value Function')
xlabel('Capital Stock')
ylabel('Lifetime Utility')

%  ... plot steady-state value
plotvdash(kstar,vstar)
plotbullet(kstar,vstar)
xlim([kmin kmax])
yl = ylim;
plottext(kstar+0.04,yl(1),'$k^*$')

% Plot shadow price function
figure
hold on
p = funeval(c,basis,k,1);
plot(k,p)
title('Shadow Price Function')
xlabel('Capital Stock')
ylabel('Shadow Price')

%  ... plot steady-state shadow price
plothdash(kstar,lstar)
plotvdash(kstar,lstar)
plotbullet(kstar,lstar)
plottext(kstar,[],'$k^*$')
plottext([],lstar,'$\lambda^*$')

% Plot residual
figure
hold on
plot(k,res)
plothdash([],0)
title('HJB Equation Residual')
xlabel('Capital Stock')
ylabel('Residual')


%% SIMULATION

% Initial state and time horizon
k0 = kmin;            % initial capital stock
T  = 50;              % time horizon

% Simulate model
[t,ksim,qsim] = docsimul(model,basis,k,q,k0,T);

% Plot simulated state and control paths
figure
hold on
plot(t,[ksim qsim])
legend('Capital Stock','Rate of Consumption')
title('Simulated Capital Stock and Rate of Consumption')
xlabel('Time')
ylabel('Quantity')

%  ... plot steady-state
plothdash([],kstar,'b')
plothdash([],qstar,'r')
plottext([],kstar,'$k^*$')
plottext([],qstar,'$q^*$')


%% SAVE FIGURES
printfigures(mfilename)


%% DOCSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, and
%   transition at an arbitrary number of ns states:
%      out = func(flag,s,x,Vs,<params>)
%   where s is ns.ds states, x is ns.dx controls, Vs is ns.ds first
%   derivatives of value function.

function out = func(flag,k,q,Vk,alpha,delta,theta)
switch flag
  case 'x'      % optimal control
    out = Vk.^(-1/theta);
  case 'f'      % reward
    out = (1/(1-theta))*q.^(1-theta);
  case 'g'      % transition
    out = k.^alpha-delta*k-q;
end