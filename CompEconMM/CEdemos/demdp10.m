%% DEMDP10 Water Resource Management Model
%
% Public authority must decide how much water to release from a reservoir so 
% as to maximize benefits derived from agricultural and recreational uses.
%
% States
%     s       reservoir level at beginning of summer
% Actions
%     x       quantity of water released for irrigation
% Parameters
%     a       producer benefit function parameters
%     b       recreational user benefit function parameters
%     delta   discount factor

function demdp10

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
a = [1 -2];                                 % producer benefit function parameters
b = [2 -3];                                 % recreational user benefit function parameters
ymean = 1.0;                                % mean rainfall
sigma = 0.2;                                % rainfall volatility
delta = 0.9;                                % discount factor

% Continuous state shock distribution
m = 3;                                            % number of shocks
[e,w] = qnwlogn(m,log(ymean)-sigma^2/2,sigma^2);  % shocks and proabilities

% Model structure
model.func = @func;                         % model functions
model.params = {a b};                       % function parameters
model.discount = delta;                     % discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 1;                               % dimension of continuous action
model.ni = 0;                               % number of discrete states
model.nj = 0;                               % number of discrete actions
model.e  = e;                               % continuous state shocks
model.w  = w;                               % continuous state shock probabilities

% Approximation structure
n    = 15;                                  % number of collocation nodes
smin =  2;                                  % minimum state
smax =  8;                                  % maximum state
basis = fundefn('cheb',n,smin,smax);        % basis functions
s = funnode(basis);                         % collocaton nodes


%% SOLUTION

% Deterministic steady-state
xstar = 1;                                  % action
sstar = 1+(a(1)*(1-delta)/b(1))^(1/b(2));   % state

% Check model derivatives
dpcheck(model,sstar,xstar)

% Compute linear-quadratic approximation at collocation nodes
[vlq,xlq] = lqapprox(model,s,sstar,xstar); 

% Solve collocation equation
[c,s,v,x,resid] = dpsolve(model,basis,vlq,xlq);


%% ANALYSIS

% Plot optimal policy
figure
plot(s,x)
title('Optimal Irrigation Policy')
xlabel('Reservoir Level')
ylabel('Irrigation')

% Plot value function
figure
plot(s,v)
title('Value Function')
xlabel('Reservoir Level')
ylabel('Social Welfare')

% Plot shadow price function
figure
pr = funeval(c,basis,s,1);
plot(s,pr)
title('Shadow Price Function')
xlabel('Reservoir Level')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
title('Bellman Equation Residual')
xlabel('Reservoir Level')
ylabel('Residual')


%% SIMULATION

% Simulation parameters
nper = 31;                              % number of periods simulated
nrep = 10000;                           % number of replications

% Initialize simulation
sinit = smin*ones(nrep,1);              % initial reservoir levels
rng('default')

% Simulate model
[ssim,xsim] = dpsimul(model,basis,nper,sinit,[],s,v,x);

% Plot simulated state path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Reservoir Level')
xlabel('Year')
ylabel('Reservoir Level')

% Plot simulated action path
figure
hold on
plot(0:nper-1,xsim(1:3,:))
plot(0:nper-1,mean(xsim),'k')
title('Simulated and Expected Irrigation')
xlabel('Year')
ylabel('Irrigation')

% Ergodic moments
ssim = ssim(:,nper);
xsim = xsim(:,nper);
savg = mean(ssim(:)); 
xavg = mean(xsim(:)); 
sstd = std(ssim(:)); 
xstd = std(xsim(:)); 
fprintf('Ergodic Moments\n') 
fprintf('                 Deterministic    Ergodic      Ergodic\n') 
fprintf('                 Steady-State      Mean     Std Deviation\n') 
fprintf('Reservoir Level     %5.2f         %5.2f         %5.2f\n'  ,[sstar savg sstd])
fprintf('Irrigation          %5.2f         %5.2f         %5.2f\n\n',[xstar xavg xstd])

% Plot ergodic state distribution
[qq,ss] = ksdensity(ssim(:),'support','positive');
figure
plot(ss,qq)
title('Ergodic Reservoir Level Distribution')
xlabel('Reservoir Level')
ylabel('Probability')
xlim([2 6])  


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

function [out1,out2,out3] = func(flag,s,x,i,j,e,a,b)

switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = s;
    out3 = [];
  case 'f'      % reward
    out1 = (a(1)/(1+a(2)))*x.^(1+a(2))+(b(1)/(1+b(2)))*(s-x).^(1+b(2));
    out2 = a(1)*x.^a(2)-b(1)*(s-x).^b(2);
    out3 = a(1)*a(2)*x.^(a(2)-1)+b(1)*b(2)*(s-x).^(b(2)-1);
  case 'g'      % transition
    out1 = s-x+e;
    out2 = -ones(size(s));
    out3 = zeros(size(s));
end