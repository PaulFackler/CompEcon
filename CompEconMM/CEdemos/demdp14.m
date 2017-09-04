%% DEMDP14 Livestock Feeding Model - Solution of Bellman Equation
%
% Farmer must decide how much to feed his livestock over a finite horizon.
% 
%  States
%      s       weight of livestock at beginning of period
%  Actions
%      x       weight gain current period
%  Parameters
%      alpha   weight gain function parameter
%      beta    weight gain function parameter
%      k       unit cost of feed
%      p       price of livestock per pound
%      N       number of feeding periods
%      s1      initial livestock weight
%      delta   discount factor

function demdp14

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.1;                           % weight gain function parameter
beta  = 2.0;                           % weight gain function parameter
k     = 0.4;                           % unit cost of feed
p     = 1.0;                           % price of livestock per pound
T     = 6;                             % number of feeding periods
s1    = 0.4;                           % initial livestock weight
delta = 0.95;                          % discount factor

% Model structure
model.horizon = T;                     % number of decision periods
model.func = @func;                    % model functions
model.discount = delta;                % discount factor
model.params = {alpha beta k};         % other parameters
model.ds = 1;                          % dimension of continuous state
model.dx = 1;                          % dimension of continuous action
model.ni = 0;                          % number of discrete states
model.nj = 0;                          % number of discrete actions

% Approximation structure
n    = 50;                             % number of collocation nodes
smax = 4.0;                            % maximum state
basis = fundefn('spli',n,s1,smax);     % basis functions
snodes = funnode(basis);               % collocaton nodes

% Intialize policy and set post-terminal value function
vterm = p*snodes;


%% SOLUTION

% Solve collocation equations
[c,s,v,x] = dpsolve(model,basis,vterm);


%% ANALYSIS

% Plot optimal policy
figure
plot(s,x)
legend('t=1','t=2','t=3','t=4','t=5','t=6')
title('Optimal Weight Gain Policy')
xlabel('Livestock Weight')
ylabel('Weight Gain')

% Plot value function
figure
plot(s,v)
legend('t=1','t=2','t=3','t=4','t=5','t=6','t=7','Location','NW')
title('Value Function')
xlabel('Livestock Weight')
ylabel('Present Value of Current and Future Profit')

% Plot shadow price function
figure
p = funeval(c,basis,s,1);
plot(s,p)
legend('t=1','t=2','t=3','t=4','t=5','t=6','t=7')
title('Shadow Price Function')
xlabel('Livestock Weight')
ylabel('Shadow Price')


%% SIMULATION

% Simulate model
[ssim,xsim] = dpsimul(model,basis,T,s1,1,s,v,x);

% Terminal weight
fprintf('Terminal Weight = %5.2f\n',ssim(end))

% Plot simulated state path
figure
plot(1:T+1,ssim)
title('Simulated Livestock Weight')
xlabel('Period')
ylabel('Weight')
xlim([1 T+1])

% Plot simulated action path
figure
plot(1:T,xsim)
title('Simulated Livestock Weight Gain')
xlabel('Period')
ylabel('Weight Gain')
xlim([1 T])


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

function [out1,out2,out3] = func(flag,s,x,~,~,~,alpha,beta,k)

ns = length(s);
switch flag
  case 'b'      % bounds
    out1 = zeros(ns,1);
    out2 = inf*ones(ns,1);
    out3 = [];
  case 'f'      % reward
    out1 = -k*(x+alpha*s).^beta;
    out2 = -k*beta*(x+alpha*s).^(beta-1);
    out3 = -k*beta*(beta-1)*(x+alpha*s).^(beta-2);
  case 'g'      % transition
    out1 = s + x;
    out2 = ones(ns,1);
    out3 = zeros(ns,1);
end