%% DEMDP08 Public Renewable Resource Model
%
% Welfare maximizing social planner must decide how much of a renewable 
% resource to harvest.
%
% States
%     s       quantity of stock available
% Actions
%     q       quantity of stock harvested
% Parameters
%     alpha   growth function parameter
%     beta    growth function parameter
%     gamma   relative risk aversion
%     kappa   unit cost of harvest
%     delta   discount factor

function demdp08

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
alpha = 4.0;                                    % growth function parameter
beta  = 1.0;                                  	% growth function parameter
gamma = 0.5;                                   	% relative risk aversion
kappa = 0.2;                                    % unit cost of harvest
delta = 0.9;                                  	% discount factor
    
% Model structure
model.func = @func;                             % model functions
model.params = {alpha beta gamma kappa};        % function parameters
model.discount = delta;                         % discount factor
model.ds = 1;                                   % dimension of continuous state
model.dx = 1;                                   % dimension of continuous action
model.ni = 0;                                   % number of discrete states
model.nj = 0;                                   % number of discrete actions

% Approximation structure
n    = 8;                                     	% number of collocation nodes
smin = 6;                                      	% minimum state
smax = 9;                                       % maximum state
basis = fundefn('cheb',n,smin,smax);            % basis functions


%% SOLUTION
  
% Steady-state
sstar = (alpha^2-1/delta^2)/(2*beta);           % stock
qstar = sstar - (delta*alpha-1)/(delta*beta); 	% action
pstar = qstar.^(-gamma);                        % market price
lstar = pstar-kappa;                            % shadow price
vstar = ((qstar.^(1-gamma))/(1-gamma)-kappa*qstar)/(1-delta); % value
fprintf('Steady States\n') 
fprintf('   Stock         %5.2f\n'  ,sstar)
fprintf('   Harvest       %5.2f\n'  ,qstar)
fprintf('   Market Price  %5.2f\n'  ,pstar)
fprintf('   Shadow Price  %5.2f\n'  ,lstar)

% Check model derivatives
dpcheck(model,sstar,qstar);

% Solve collocation equation
[c,s,v,q,resid] = dpsolve(model,basis);


%% ANALYSIS

% Plot optimal policy
figure
hold on
plot(s,q)
title('Optimal Harvest Policy')
xlabel('Stock')
ylabel('Quantity Harvested')

% ...plot steady-state action
plotbullet(sstar,qstar)
plotvdash(sstar,qstar)
plothdash(sstar,qstar)
plottext(sstar+0.02,[],'$s^*$')
plottext([],qstar,'$q^*$')

% Plot value function
figure
hold on
plot(s,v)
title('Value Function')
xlabel('Stock')
ylabel('Social Welfare')

%  ... plot steady-state value
plotvdash(sstar,vstar)
plotbullet(sstar,vstar)
plottext(sstar+0.01,[],'$s^*$')

% Plot shadow price function
figure
hold on
p = funeval(c,basis,s,1);
plot(s,p)
title('Shadow Price Function')
xlabel('Stock')
ylabel('Shadow Price')

% ...plot steady-state shadow price
plotvdash(sstar,lstar)
plothdash(sstar,lstar)
plotbullet(sstar,lstar)
plottext(sstar+0.02,[],'$s^*$')
plottext([],lstar,'$\lambda^*$')

% Plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
title('Bellman Equation Residual')
xlabel('Stock')
ylabel('Residual')


%% SIMULATION

% Simulation parameters
nper = 11;                              % number of periods simulated

% Initialize simulation
sinit = smin;                           % agent possesses minimal wealth

% Simulate model
[ssim,qsim] = dpsimul(model,basis,nper,sinit,[],s,v,q);

% Plot simulated state and policy paths
figure
hold on
plot(0:nper-1,ssim,0:nper-1,qsim)
legend('Stock','Harvest')
title('Simulated Stock and Harvest')
xlabel('Period')
ylabel('Quantity')

% ...plot steady-state
plothdash([],sstar,'b')
plothdash([],qstar,'r')
plottext([],sstar,'$s^*$')
plottext([],qstar,'$q^*$')


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

function [out1,out2,out3] = func(flag,s,q,~,~,~,alpha,beta,gamma,kappa)

switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = s;
    out3 = [];
  case 'f'      % reward
    out1 = (q.^(1-gamma))/(1-gamma)-kappa*q;
    out2 = q.^(-gamma)-kappa;
    out3 = -gamma*q.^(-gamma-1);
  case 'g'      % transition
    out1 = alpha*(s-q) - 0.5*beta*(s-q).^2;
    out2 = -alpha + beta*(s-q);
    out3 = -beta*ones(size(s));
end