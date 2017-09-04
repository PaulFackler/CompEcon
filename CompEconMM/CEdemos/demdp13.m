%% DEMDP13 Inventory Management Model
%
% Profit maximizing entrepeneur must decide how much to produce and how 
% much inventory to hold.
%
% States
%     s1      market price
%     s2      initial stock of inventory
% Actions
%     x1      quantity produced
%     x2      initial stock of inventory
% Parameters
%     c       production cost function parameters
%     k       inventory holding cost function parameters
%     pbar	  long-run mean price
%     rho     mean-reversion coefficient
%     sigma   standard deviation of price shocks
%     delta   discount factor

function demdp13

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
c     = [0.5 0.1];                      % production cost function parameters
k     = [0.1 0.1];                      % inventory holding cost function parameters
pbar  = 1.0;                            % long-run mean price
rho   = 0.5;                            % mean-reversion coefficient
sigma = 0.2;                            % standard deviation of price shocks
delta = 0.9;                            % discount factor

% Continuous state shock distribution
m   = 3;                                % number of shocks
[e,w] = qnwnorm(m,0,sigma^2);           % shocks and probabilities

% Model structure
model.func = @func;                     % model functions
model.params = {c k pbar rho};          % function parameters
model.discount = delta;                 % discount factor
model.ds = 2;                           % dimension of continuous state
model.dx = 2;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities

% Approximation structure
n = [15 60];                           	% number of collocation node coordinates, per dimension
smin = [0.2 0.0];                     	% minimum states
smax = [2.5 3.0];                     	% maximum states
basis = fundefn('spli',n,smin,smax);    % basis functions


%% SOLUTION

% Deterministic steady-state
sstar = [pbar 0];                       % state
xstar = [(pbar-c(1))/c(2) 0];           % action

% Check model derivatives
dpcheck(model,sstar,xstar)

% Solve collocation equation
[c,s,v,x,resid,n] = dpsolve(model,basis);


%% ANALYSIS

% Reshape output for plotting
s1 = reshape(s(:,1),n);
s2 = reshape(s(:,2),n);
v = reshape(v,n);
x1 = reshape(x(:,1),n);
x2 = reshape(x(:,2),n);
resid = reshape(resid,n);

% Shadow prices
p1 = funeval(c,basis,s,[1 0]);
p1 = reshape(p1,n);

% Plot optimal policy 1
figure
mesh(s1,s2,x1)
title('Optimal Production Policy')
xlabel('Market Price')
ylabel('Initial Inventory')
zlabel('Production')

% Plot optimal policy 2
figure
mesh(s1,s2,x2)
title('Optimal Inventory Policy')
xlabel('Market Price')
ylabel('Initial Inventory')
zlabel('Ending Inventory')

% Plot value function
figure
mesh(s1,s2,v)
title('Value Function')
xlabel('Market Price')
ylabel('Initial Inventory')
zlabel('Value of the Firm')

% Plot shadow price function 1
figure
mesh(s1,s2,p1)
title('Shadow Price of Inventories')
xlabel('Market Price')
ylabel('Initial Inventory')
zlabel('Price')

% Plot residual
figure
mesh(s1,s2,resid)
title('Bellman Equation Residual')
xlabel('Market Price')
ylabel('Initial Inventory')
zlabel('Residual')


%% SIMULATION

% Simulation parameters
nper = 26;                              % number of periods simulated
nrep = 50000;                           % number of replications

% Initialize simulation
sinit = [pbar*ones(nrep,1) zeros(nrep,1)];
rng('default')

% Simulate model
[ssim,xsim] = dpsimul(model,basis,nper,sinit,[],s,v,x);
s1sim = ssim(:,:,1);
s2sim = ssim(:,:,2);
x1sim = xsim(:,:,1);

% Ergodic moments
s1avg = mean(s1sim(:)); 
s2avg = mean(s2sim(:)); 
x1avg = mean(x1sim(:)); 
s1std = std(s1sim(:)); 
s2std = std(s2sim(:)); 
x1std = std(x1sim(:)); 
fprintf('Ergodic Moments\n') 
fprintf('              Deterministic    Ergodic      Ergodic\n') 
fprintf('              Steady-State      Mean     Std Deviation\n') 
fprintf('Market Price     %5.3f         %5.3f         %5.3f\n',[sstar(1) s1avg s1std])
fprintf('Inventory        %5.3f         %5.3f         %5.3f\n',[sstar(2) s2avg s2std])
fprintf('Production       %5.3f         %5.3f         %5.3f\n',[xstar(1) x1avg x1std])

% Plot simulated and expected state paths 1
figure
hold on
plot(0:nper-1,s1sim(1:3,:))
plot(0:nper-1,mean(s1sim),'k')
plot(nper-1,s1avg,'k*')
title('Simulated and Expected Market Price')
xlabel('Period')
ylabel('Market Price')

% Plot simulated and expected state paths 2
figure
hold on
plot(0:nper-1,s2sim(1:3,:))
plot(0:nper-1,mean(s2sim),'k')
plot(nper-1,s2avg,'k*')
title('Simulated and Expected Ending Inventory')
xlabel('Period')
ylabel('Ending Inventory')

% Plot simulated and expected action paths
figure
hold on
plot(0:nper-1,squeeze(x1sim(1:3,:)))
plot(0:nper-1,mean(x1sim),'k')
plot(nper-1,x1avg,'k*')
title('Simulated and Expected Production')
xlabel('Period')
ylabel('Production')


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

function [out1,out2,out3] = func(flag,s,x,~,~,e,c,k,pbar,rho)

n = size(s,1);
ds = 2;
dx = 2;

switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = inf*ones(size(s));
    out3 = [];
  case 'f'      % reward
    out2 = zeros(n,dx);
    out3 = zeros(n,dx,dx);
    out1 = s(:,1).*(s(:,2)+x(:,1)-x(:,2)) ...
      - (c(1)+0.5*c(2)*x(:,1)).*x(:,1) ...
      - (k(1)+0.5*k(2)*x(:,2)).*x(:,2);
    out2(:,1) =  s(:,1) - (c(1)+c(2)*x(:,1));
    out2(:,2) = -s(:,1) - (k(1)+k(2)*x(:,2));
    out3(:,1,1) = -c(2)*ones(n,1);
    out3(:,2,2) = -k(2)*ones(n,1);
  case 'g'      % transition
    out2 = zeros(n,ds,dx);
    out3 = zeros(n,ds,dx,dx);
    out1 = [pbar+rho*(s(:,1)-pbar)+e x(:,2)];
    out2(:,2,2) = ones(n,1);
end