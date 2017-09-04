%% DEMDP06SHORT  Deterministic Ramsey Optimal Economic Growth Model
%
% Welfare maximizing social planner must decide how much society should
% consume and invest.  Model is of special interest because it has known 
% closed-form solution.
% 
% States
%     s       stock of wealth
% Actions
%     k       capital investment
% Parameters
%     beta	  capital production elasticity
%     delta   discount factor

function demdp06short

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
beta  = 0.4;                            % capital production elasticity
delta = 0.9;                            % discount factor

% Model structure
model.func = @func;                     % model functions
model.params = {beta};                  % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions

% Approximation structure
n     = 15;                             % number of collocation nodes
smin  = 0.2;                            % minimum state
smax  = 1.0;                            % maximum state
basis = fundefn('cheb',n,smin,smax);    % basis functions
s = funnode(basis);                     % collocaton nodes


%% SOLUTION

% Steady-state
sstar = (beta*delta)^(beta/(1-beta));   % wealth
kstar = beta*delta*sstar;               % capital investment
vstar = log(sstar-kstar)/(1-delta);     % value
pstar = 1/(sstar*(1-beta*delta));       % shadow price
fprintf('Steady States\n') 
fprintf('   Wealth                %5.2f\n'  ,sstar)
fprintf('   Capital Investment    %5.2f\n'  ,kstar)
fprintf('   Shadow Price          %5.2f\n'  ,pstar)
  
% Analytic solution at collocation nodes
b = 1/(1-delta*beta);
vtrue = vstar + b*(log(s)-log(sstar));
ktrue = delta*beta*s;

% Solve collocation equation
[c,s,v,k,resid] = dpsolve(model,basis,vtrue,ktrue);
   
% Linear-Quadratic approximation at refined state grid
[vlq,klq,plq] = lqapprox(model,s,sstar,kstar);

% Analytic solution at refined state grid
vtrue = vstar + b*(log(s)-log(sstar));


%% ANALYSIS

% Plot optimal policy
figure
hold on
plot(s,[k klq])
legend('Chebychev Collocation','L-Q Approximation')
title('Optimal Investment Policy')
xlabel('Wealth')
ylabel('Capital Investment')

%  ... plot steady-state action
plothdash(sstar,kstar)
plotvdash(sstar,kstar)
plotbullet(sstar,kstar)
plottext(sstar+0.01,[],'$s^*$')
plottext([],kstar,'$k^*$')

% Plot value function
figure
hold on
plot(s,[v vlq])
legend('Chebychev Collocation','L-Q Approximation')
title('Value Function')
xlabel('Wealth')
ylabel('Lifetime Utility')

%  ... plot steady-state value
plotvdash(sstar,vstar)
plotbullet(sstar,vstar)
plottext(sstar+0.01,[],'$s^*$')

% Plot shadow price function
figure
hold on
pr = funeval(c,basis,s,1);
plot(s,[pr plq])
legend('Chebychev Collocation','L-Q Approximation')
title('Shadow Price Function')
xlabel('Wealth')
ylabel('Shadow Price')

%  ... plot steady-state shadow price
plothdash(sstar,pstar)
plotvdash(sstar,pstar)
plotbullet(sstar,pstar)
plottext(sstar+0.01,[],'$s^*$')
plottext([],pstar-0.6,'$\lambda^*$')

% Plot Chebychev collocation Bellman equation residual and approximation error
figure
hold on
plot(s,[resid v-vtrue])
plothdash([],0)
legend('Residual','Error')
title('Chebychev Collocation Residual and Approximation Error')
xlabel('Wealth')
ylabel('Residual/Error')

% Plot Linear-Quadratic approximation error
figure
plot(s,vlq-vtrue)
title('Linear-Quadratic Approximation Error')
xlabel('Wealth')
ylabel('Error')


%% SIMULATION

% Simulation parameters
nper = 21;                              % number of periods simulated
nrep = 1;                               % number of replications

% Initialize simulation
sinit = smax*ones(nrep,1);              % initial wealths
rng('default')

% Simulate model
[ssim,ksim] = dpsimul(model,basis,nper,sinit,[],s,v,k);

% Plot simulated state and action paths
figure
hold on
plot(0:nper-1,ssim,0:nper-1,ksim)
legend('Wealth','Investment')
title('Simulated Wealth and Investment')
xlabel('Period')
ylabel('Amount')

%  ... plot steady-state
plothdash([],sstar,'b')
plothdash([],kstar,'r')
plottext([],sstar,'$s^*$')
plottext([],kstar,'$k^*$')


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

function [out1,out2,out3] = func(flag,s,k,~,~,~,beta)

n = length(s);
switch flag
  case 'b'      % bounds
    out1 = zeros(n,1);
    out2 = s;
    out3 = [];
  case 'f'      % reward
    out1 = log(s-k);
    out2 = -(s-k).^(-1);
    out3 = -(s-k).^(-2);
  case 'g'      % transition
    out1 = k.^beta;
    out2 = beta*k.^(beta-1);
    out3 = (beta-1)*beta*k.^(beta-2);
end