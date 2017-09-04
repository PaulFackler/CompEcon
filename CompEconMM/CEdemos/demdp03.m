%% DEMDP03 Dixit Industry Entry-Exit Model
%
% Profit maximizing firm must decide whether to operate or shut down, given
% its current operational state and current profitability, subject to
% transactions costs.
%
% States
%     p       current profitability
%     i       active (1) or idle (2) last period
% Actions
%     j       active (1) or idle (2) this period
% Parameters
%     pbar    long-run mean profitability
%     gamma   profitability autoregressive coefficient
%     kappa   cost of reopenning idle firm
%     sigma   standard deviation of profitability shock
%     delta   discount factor

function demdp03
  
% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
pbar  = 1.0;                            % long-run mean profitability
gamma = 0.7;                            % profitability autoregressive coefficient
kappa =  10;                            % cost of reopenning idle firm
sigma = 1.0;                            % standard deviation of profitability shock
delta = 0.9;                            % discount factor
  
% Continuous state shock distribution
m = 5;                                  % number of shocks
[e,w] = qnwnorm(m,0,sigma.^2);          % shocks and probabilities
  
% Model structure
model.func = @func;                     % model functions
model.params = {pbar gamma kappa};      % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 0;                           % dimension of continuous action
model.ni = 2;                           % number of discrete states
model.nj = 2;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities
model.h  = [1 1;2 2];                   % deterministic discrete state transitions

% Approximation structure
n = 250;                                % number of collocation nodes
pmin = -20;                             % minimum state
pmax =  20;                             % maximum state
basis = fundefn('spli',n,pmin,pmax);    % basis functions


%% SOLUTION

% Solve collocation equation
[~,sr,vr,~,resid] = dpsolve(model,basis); 

% Critical profitabilities and values
pcrit = [interp1(vr(:,1,1)-vr(:,1,2),sr,0); interp1(vr(:,2,1)-vr(:,2,2),sr,0)];
vcrit = [interp1(sr,vr(:,1,1),pcrit(1)); interp1(sr,vr(:,2,1),pcrit(2))];
fprintf('Critical Profitabilities\n')
fprintf('   Exit  %5.2f\n'    ,pcrit(1))  
fprintf('   Entry %5.2f\n\n'  ,pcrit(2))


%% ANALYSIS

% Plot action-contingent value functions
figure
hold on
plot(sr,[vr(:,1,1) vr(:,2,1) vr(:,2,2)])
legend('Keep Active Firm Open','Reopen Idle Firm','Shut Down')
title('Action-Contingent Value Functions')
xlabel('Profitability')
ylabel('Value of Firm')

%  ... plot critical profitabilities
for i=1:2
  plotvdash(pcrit(i),vcrit(i))
  plotbullet(pcrit(i),vcrit(i))
  plottext(pcrit(i)+0.6,-40,['$p^*_' int2str(2-i) '$'])
end

% Plot residual
figure
hold on
plot(sr,100*resid./max(vr,[],3))
plothdash([],0)
legend('Active','Idle')
title('Bellman Equation Residual')
xlabel('Profitability')
ylabel('Percent Residual')

%  ... plot critical profitabilities
for i=1:2
  plotvdash(pcrit(i),0)
  plotbullet(pcrit(i),0)
  plottext(pcrit(i)+0.6,[],['$p^*_' int2str(2-i) '$'])
end


%% SIMULATION

% Simulation parameters
nper = 51;                              % number of periods simulated
nrep = 50000;                           % number of replications

% Initialize simulation
pinit = pbar*ones(nrep,1);              % initial profitabilities
iinit = ones(nrep,1);                   % initial firm activity states
rng('default')

% Simulate Model
[ssim,~,isim] = dpsimul(model,basis,nper,pinit,iinit,sr,vr);

% Ergodic moments
isim = 2 - isim; % Convert to 0=idle, 1=active
savg = mean(ssim(:,nper));
iavg = mean(isim(:,nper));
sstd = std(ssim(:,nper));
istd = std(isim(:,nper));
fprintf('Ergodic Moments\n') 
fprintf('                Ergodic      Ergodic\n') 
fprintf('                 Mean     Std Deviation\n') 
fprintf('Profitability   %5.3f         %5.3f\n'  ,[savg sstd])
fprintf('Activity        %5.3f         %5.3f\n\n',[iavg istd])

% Plot simulated and expected continuous state path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Profitabilities')
xlabel('Period')
ylabel('Profitability')

% Plot expected discrete state path
figure
hold on
plot(0:nper-1,mean(isim),'k')
title('Probability of Operation')
xlabel('Period')
ylabel('Probability')


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

function out = func(flag,p,~,i,j,e,pbar,gamma,kappa)

switch flag
  case 'f'      % reward
    out = p.*(j==1) - kappa.*(i==2).*(j==1);
  case 'g'      % transition
    out = pbar+gamma*(p-pbar)+e;
end