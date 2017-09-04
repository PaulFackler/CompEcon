%% DEMDP01C Timber Harvesting Model - Cubic Spline Approximation Using DPSOLVE
%
% Profit maximizing a commercial tree stand owner must decide when to
% clear-cut and replant.  This program uses a cubic spline approximation
% for the value function and solves the collocation equation direclly using
% the CompEcon Toolbox routine DPSOLVE.
%
% States
%     s       stand biomass
% Actions
%     j       clear cut/replant (2), not clear cut (1)
% Parameters
%     price   unit price of biomass
%     kappa   clearcut-replant cost
%     smax    stand carrying capacity
%     gamma   biomass growth parameter
%     delta   discount factor

function demdp01c

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
price = 1.0;                           	% unit price of biomass
kappa = 0.2;                            % clearcut-replant cost
smax  = 0.5;                            % stand carrying capacity
gamma = 0.1;                            % biomass growth parameter
delta = 0.9;                            % discount factor

% Growth function
h = @(s) s+gamma*(smax-s);
  
% Model structure
model.func = @func;                     % model functions
model.params = {price kappa h};         % function parameters
model.discount = delta;                	% discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 0;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 2;                           % number of discrete actions

% Approximation structure
n = 200;                                % number of collocation nodes
basis = fundefn('spli',n,0,smax);       % basis functions


%% SOLUTION

% Solve collocation equation
[~,sr,vr,~,resid] = dpsolve(model,basis); 

% Critical biomass and value
scrit = interp1(vr(:,1)-vr(:,2),sr,0);
vcrit = interp1(sr,vr(:,1),scrit);
fprintf('Critical Biomass           %9.3f\n',scrit) 


%% ANALYSIS

% Plot action-contingent value functions
figure
hold on
plot(sr,vr)
legend('Grow','Clear-Cut')
title('Action-Contingent Value Functions')
xlabel('Biomass')
ylabel('Value of Stand')

% ... plot critical biomass
plotvdash(scrit,vcrit)
plotbullet(scrit,vcrit)
plottext(scrit+0.01,[],'$s^*$')

% Plot residual
figure
hold on
plot(sr,100*resid./max(vr,[],2))
plothdash([],0)
title('Bellman Equation Residual')
xlabel('Biomass')
ylabel('Percent Residual')

% ... plot critical biomass
plotvdash(scrit,0)
plotbullet(scrit,0)
plottext(scrit+0.01,[],'$s^*$')


%% SIMULATION

% Simulation parameters
nper = 31;                              % number of periods simulated
time = 0:nper-1;                        % periods simulated

% Initialize simulation
s = 0;                                  % initial biomass

% Simulate model
ssim = zeros(nper,1);
for ip=1:nper
  ssim(ip) = s;
  if s<scrit,
    s = h(s);
  else
    s = 0;
  end
end

% Compute mean annual harvest and optimal rotation cycle
s = 0;
for n=1:100
  if s>scrit, break, end
  s = h(s);
end
fprintf('Mean Annual Harvest        %9.3f\n',s/n) 
fprintf('Rotation Cycle in Years    %9i\n\n',n) 

% Plot state path
figure
plot(time,ssim)
title('Simulated Biomass')
xlabel('Period')
ylabel('Biomass')


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

function out = func(flag,s,~,~,j,~,price,kappa,h)
n = length(s);
switch flag
  case 'f'      % reward
    out = (price*s-kappa).*(j-1);
  case 'g'      % transition
    if j==1
      out = h(s);
    else
      out = h(0)*ones(n,1);
    end
end