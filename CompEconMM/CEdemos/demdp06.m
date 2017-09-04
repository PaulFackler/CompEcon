%% DEMDP06  Stochastic Ramsey Optimal Economic Growth Model
%
% Welfare maximizing social planner must decide how much society should
% consume and invest.  Model is of special interest because it has a known 
% closed-form solution.
% 
% States
%     s       stock of wealth
% Actions
%     k       capital investment
% Parameters
%     beta	  capital production elasticity
%     delta   discount factor
%     sigma   production shock volatility

function demdp06

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
beta  = 0.5;                            % capital production elasticity
delta = 0.9;                            % discount factor
sigma = 0.0;                            % production shock volatility

% Continuous state shock distribution
if sigma>0
  m = 7;                                  % number of shocks
  [e,w] = qnwlogn(m,-sigma^2/2,sigma^2); 	% shocks and probabilities
else
  e = 1;
  w = 1;
end

% Model structure
model.func = @func;                     % model functions
model.params = {beta};                  % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities

% Approximation structure
n     = 15;                             % number of collocation nodes
smin  = 0.2;                            % minimum state
smax  = 1.0;                            % maximum state
basis = fundefn('cheb',n,smin,smax);    % basis functions
s = funnode(basis);                     % collocaton nodes


%% SOLUTION

% Deterministic steady-state
sstar = (beta*delta)^(beta/(1-beta));   % wealth
kstar = beta*delta*sstar;               % capital investment
vstar = log(sstar-kstar)/(1-delta);     % value
pstar = 1/(sstar*(1-beta*delta));       % shadow price

% Analytic solution at collocation nodes
b = 1/(1-delta*beta);
vtrue = vstar - (0.5*b*delta*sigma^2)/(1-delta) + b*(log(s)-log(sstar));
ktrue = delta*beta*s;

% Solve collocation equation
[c,s,v,k,resid] = dpsolve(model,basis,vtrue,ktrue);
   
% Linear-Quadratic approximation at refined state grid
[vlq,klq,plq] = lqapprox(model,s,sstar,kstar);

% Analytic solution at refined state grid
vtrue = vstar - (0.5*b*delta*sigma^2)/(1-delta) + b*(log(s)-log(sstar));


%% ANALYSIS

% Plot optimal policy
figure
hold on
plot(s,[k klq])
legend('Chebychev Collocation','L-Q Approximation')
title('Optimal Investment Policy')
xlabel('Wealth')
ylabel('Capital Investment')

%  ... plot deterministic steady-state action
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

%  ... plot deterministic steady-state value
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

%  ... plot deterministic steady-state shadow price
plothdash(sstar,pstar)
plotvdash(sstar,pstar)
plotbullet(sstar,pstar)
plottext(sstar+0.01,[],'$s^*$')
plottext([],pstar-1.0,'$\lambda^*$')

% Plot Chebychev collocation Bellman equation residual and approximation error
figure
hold on
plot(s,[resid v-vtrue])
plothdash([],0)
title('Chebychev Collocation Residual and Approximation Error')
legend('Residual','Error')
xlabel('Wealth')
ylabel('Residual/Error')

% Plot Linear-Quadratic approximation error
figure
hold on
plot(s,vlq-vtrue)
title('Linear-Quadratic Approximation Error')
xlabel('Wealth')
ylabel('Error')


%% SIMULATION

% Simulation parameters
nper = 21;                              % number of periods simulated
if sigma>0
  nrep = 10000;                         % number of replications
else
  nrep = 1;                             % number of replications
end

% Initialize simulation
sinit = smin*ones(nrep,1);              % initial wealths
rng('default')

% Simulate model
[ssim,ksim] = dpsimul(model,basis,nper,sinit,[],s,v,k);

% Plot simulated state paths
figure
hold on
if sigma>0
  plot(0:nper-1,ssim(1:3,:))
  plot(0:nper-1,mean(ssim),'k')
else
  plot(0:nper-1,ssim)
end
title('Simulated and Expected Wealth')
xlabel('Period')
ylabel('Wealth')

%  ... plot steady-state
plothdash([],sstar)
plottext([],sstar,'$s^*$')

% Plot simulated action paths
figure
hold on
if sigma>0
  plot(0:nper-1,ksim(1:3,:))
  plot(0:nper-1,mean(ksim),'k')
else
  plot(0:nper-1,ksim)
end
title('Simulated and Expected Investment')
xlabel('Period')
ylabel('Investment')

%  ... plot steady-state
plothdash([],kstar)
plottext([],kstar,'$k^*$')

% Ergodic moments
ssim = ssim(:,nper);
ksim = ksim(:,nper);
savg = mean(ssim(:));
kavg = mean(ksim(:));
sstd = std(ssim(:));
kstd = std(ksim(:));
fprintf('Ergodic Moments\n') 
fprintf('          Deterministic    Ergodic      Ergodic\n') 
fprintf('          Steady-State      Mean     Std Deviation\n') 
fprintf('Wealth       %5.3f         %5.3f         %5.3f\n'  ,[sstar savg sstd])
fprintf('Investment   %5.3f         %5.3f         %5.3f\n\n',[kstar kavg kstd])

% Plot ergodic state distribution
if sigma>0
  [qq,ss] = ksdensity(ssim(:),'support','positive');
  figure
  plot(ss,qq)
  title('Ergodic Wealth Distribution')
  xlabel('Wealth')
  ylabel('Probability')
  xlim([0.2 0.8])
end


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

function [out1,out2,out3] = func(flag,s,k,~,~,e,beta)

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
    out1 = e.*k.^beta;
    out2 = beta*e.*k.^(beta-1);
    out3 = (beta-1)*beta*e.*k.^(beta-2);
end