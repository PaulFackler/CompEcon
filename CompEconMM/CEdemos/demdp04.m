%% DEMDP04 Job Search Model
%
% Infinitely-lived worker must decide whether to quit, if employed, or 
% search for a job, if unemployed, given prevailing market wages.
%
% States
%     w       prevailing wage
%     i       unemployed (1) or employed (2) at beginning of period
% Actions
%     j       idle (1) or active (i.e., work or search) (2) this period
% Parameters
%     v        benefit of pure leisure
%     wbar     long-run mean wage
%     gamma    wage reversion rate
%     p0       probability of finding job
%     p1       probability of keeping job
%     sigma    standard deviation of wage shock
%     delta    discount factor

function demdp04

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
u     =  90;                            % unemployment benefit
v     =  95;                            % benefit of pure leisure
wbar  = 100;                           	% long-run mean wage
gamma = 0.40;                           % wage reversion rate
p0    = 0.20;                         	% probability of finding job
p1    = 0.90;                           % probability of keeping job
sigma = 5;                             	% standard deviation of wage shock
delta = 0.95;                           % discount factor
   
% Continuous state shock distribution
m = 15;                                	% number of shocks
[e,w] = qnwnorm(m,0,sigma^2);           % shocks and probabilities

% Stochastic discrete state transition probabilities
q = zeros(2,2,2);
q(1,2,2) = p0;
q(2,2,2) = p1;
q(:,1,:) = 1-q(:,2,:);

% Model structure
model.func = @func;                     % model functions
model.params = {u v wbar gamma};      	% function parameters
model.discount = delta;                	% discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 0;                           % dimension of continuous action
model.ni = 2;                           % number of discrete states
model.nj = 2;                           % number of discrete actions
model.e  = e;                          	% continuous state shocks
model.w  = w;                          	% continuous state shock probabilities
model.q  = q;                          	% discrete state transition probabilities

% Approximation structure
n = 150;                                % number of collocation nodes
wmin =   0;                             % minimum state
wmax = 200;                             % maximum state
basis = fundefn('spli',n,wmin,wmax);    % basis functions


%% SOLUTION

% Solve collocation equation
[~,sr,vr,~,resid] = dpsolve(model,basis); 

% Critical action wages
wcrit1 = interp1(vr(:,1,1)-vr(:,1,2),sr,0);
vcrit1 = interp1(sr,vr(:,1,1),wcrit1);
wcrit2 = interp1(vr(:,2,1)-vr(:,2,2),sr,0);
vcrit2 = interp1(sr,vr(:,2,1),wcrit2);
fprintf('Critical Wages\n') 
fprintf('   Search  %5.1f\n',wcrit1) 
fprintf('   Quit    %5.1f\n',wcrit2)


%% ANALYSIS

% Plot action-contingent value function - unemployed
figure
hold on
plot(sr,squeeze(vr(:,1,:)))
legend('Do Not Search','Search')
title('Action-Contingent Value Function - Unemployed')
xlabel('Wage')
ylabel('Lifetime Utility')

%  ... plot critical search wage
plotvdash(wcrit1,vcrit1)
plotbullet(wcrit1,vcrit1)
plottext(wcrit1+1,[],'$w_0^*$')

% Plot action-contingent value function - employed
figure
hold on
plot(sr,squeeze(vr(:,2,:)))
legend('Quit','Work')
title('Action-Contingent Value Function - Employed')
xlabel('Wage')
ylabel('Lifetime Utility')

%  ... plot critical quit wage
plotvdash(wcrit2,vcrit2)
plotbullet(wcrit2,vcrit2)
plottext(wcrit2+1,[],'$w_1^*$')

% Plot residual
figure
hold on
plot(sr,100*resid./max(vr,[],3))
plothdash([],0)
legend('Unemployed','Employed')
title('Bellman Equation Residual')
xlabel('Wage')
ylabel('Percent Residual')

%  ... plot critical wages
plotvdash(wcrit1,0)
plotbullet(wcrit1,0)
plottext(wcrit1+1,[],'$w_0^*$')
plotvdash(wcrit2,0)
plotbullet(wcrit2,0)
plottext(wcrit2+1,[],'$w_1^*$')


%% SIMULATION

% Simulation parameters
nper = 41;                              % number of periods simulated
nrep = 10000;                           % number of replications

% Initialize simulation
sinit = wbar*ones(nrep,1);              % initial wages
iinit = ones(nrep,1);                   % initial employment states
rng('default')

% Simulate model
[ssim,~,isim] = dpsimul(model,basis,nper,sinit,iinit,sr,vr);

% Ergodic moments
isim = isim-1;
savg = mean(ssim(:,nper));
iavg = mean(isim(:,nper));
sstd = std(ssim(:,nper));
istd = std(isim(:,nper));
fprintf('\nErgodic Moments\n') 
fprintf('                Ergodic      Ergodic\n') 
fprintf('                 Mean     Std Deviation\n') 
fprintf('Wage           %6.2f         %6.2f\n',[savg sstd])
fprintf('Unemployment   %6.2f         %6.2f\n',[1-iavg istd])

% Plot expected discrete state path
figure
hold on
plot(0:nper-1,1-mean(isim),'k')
title('Probability of Unemployment')
xlabel('Period')
ylabel('Probability')

% Plot simulated and expected continuous state path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Wage')
xlabel('Period')
ylabel('Wage')


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

function out = func(flag,w,~,i,j,e,u,v,wbar,gamma)

switch flag
  case 'f'      % reward
    out = (j==1)*v + (j==2)*((i==1)*u+(i==2)*w);
  case 'g'      % transition
    out = wbar+gamma*(w-wbar)+e;
end