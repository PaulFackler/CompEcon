%% DEMDP02 Asset Replacement Model
%
% Profit-maximizing manufacturer must decide when to replace an aging asset.
%
% States
%     p       unit profit contribution
%     a       asset age (1..A)
% Actions
%     j       keep(1) or replace(2) asset
% Parameters
%     A       maximum asset age 
%     alpha   production function coefficients
%     kappa   net replacement cost
%     pbar    long-run mean unit profit contribution
%     gamma   unit profit contribution autoregression coefficient
%     sigma   standard deviation of unit profit contribution shock
%     delta   discount factor

function demdp02

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
A       = 6;                                % maximum asset age 
alpha   = [50 -2.5 -2.5];                   % production function coefficients
kappa   = 40;                               % net replacement cost
pbar    = 1;                                % long-run mean unit profit contribution
gamma   = 0.5;                              % unit profit contribution autoregression coefficient
sigma   = 0.15;                             % standard deviation of unit profit contribution shock
delta   = 0.9;                              % discount factor 

% Continuous state shock distribution
m = 5;                                      % number of shocks
[e,w] = qnwnorm(m,0,sigma^2);               % shocks and probabilities

% Deterministic discrete state transitions
h = [2:A 1; ones(1,A)];

% Model structure
model.func = @func;                         % model functions
model.params = {A alpha kappa pbar gamma};  % function parameters
model.discount = delta;                   	% discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 0;                               % dimension of continuous action
model.ni = A;                               % number of discrete states
model.nj = 2;                               % number of discrete actions
model.e  = e;                              	% continuous state shocks
model.w  = w;                              	% continuous state shock probabilities
model.h  = h;                              	% deterministic discrete state transitions

% Approximation structure
n  = 200;                                   % number of collocation nodes
pmin = 0;                                   % minimum state
pmax = 2;                                   % maximum state
basis = fundefn('spli',n,pmin,pmax);        % basis functions


%% SOLUTION

% Solve collocation equation
[~,pr,vr,~,resid] = dpsolve(model,basis); 

% Critical unit profit contributions and values
fprintf('Critical Replacement Unit Profit Contributions\n')
pcrit = zeros(A,1);
vcrit = zeros(A,1);
for a=1:A-1
  pcrit(a) = interp1(vr(:,a,1)-vr(:,a,2),pr,0);
  vcrit(a) = interp1(pr,vr(:,a,1),pcrit(a));
  if isnan(pcrit), continue, end
  fprintf('   Age %2i  %5.2f\n'  ,a,pcrit(a))
end


%% ANALYSIS

% Plot action-contingent value functions
figure
hold on
plot(pr,squeeze(vr(:,:,1)))
legend('Keep a=1','Keep a=2','Keep a=3','Keep a=4','Keep a=5','Replace')
title('Action-Contingent Value Functions')
xlabel('Unit Profit Contribution')
ylabel('Value of Production Operation')

% ... plot critical unit profit contributions
for a=1:A-1
  plotvdash(pcrit(a),vcrit(a))
  plotbullet(pcrit(a),vcrit(a))
  plottext(pcrit(a)+0.01,[],['$p^*_' int2str(a) '$'])
end

% Plot residual
figure
hold on
plot(pr,100*resid./max(vr,[],3))
plothdash([],0)
legend('a=1','a=2','a=3','a=4','a=5','a=6','Location','NE')
title('Bellman Equation Residual')
xlabel('Unit Profit Contribution')
ylabel('Percent Residual')

% ... plot critical unit profit contributions
for a=1:A-1
  plotvdash(pcrit(a),0)
  plotbullet(pcrit(a),0)
  plottext(pcrit(a)+0.01,[],['$p^*_' int2str(a) '$'])
end


%% SIMULATION

% Simulation parameters
nper = 51;                              % number of periods simulated
nrep = 10000;                           % number of replications

% Initialize simulation
sinit = pbar*ones(nrep,1);              % initial unit profit contributions
iinit = ones(nrep,1);                   % initial asset ages
rng('default')

% Simulate model
[ssim,~,isim,~] = dpsimul(model,basis,nper,sinit,iinit,pr,vr);

% Ergodic moments
savg = mean(ssim(:));
iavg = mean(isim(:));
sstd = std(ssim(:));
istd = std(isim(:));
fprintf('\nErgodic Moments\n') 
fprintf('                           Mean     Std Deviation\n') 
fprintf('Unit Profit Contribution  %5.3f         %5.3f\n'  ,[savg sstd])
fprintf('Age                       %5.3f         %5.3f\n\n',[iavg istd])

% Plot simulated and expected continuous state path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Unit Profit Contribution')
xlabel('Period')
xlabel('Unit Profit Contribution')


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

function out = func(flag,p,~,a,j,e,A,alpha,kappa,pbar,gamma)

switch flag
  case 'f'      % reward
    if j==2||a==A
      out = p*50-kappa;
    else
      out = p*(alpha(1)+alpha(2)*a+alpha(3)*a.^2);
    end
  case 'g'      % transition
    out = pbar + gamma*(p-pbar) + e;
end