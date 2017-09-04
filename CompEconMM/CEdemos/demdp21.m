%% DEMDP21 Savings and Insurance
%
% An infinitely-lived agent begins each period endowed with predetermined
% wealth w, which he must allocate among consumption, savings, and
% purchases of insurance. The agent receives a certain income ybar the
% following period, but can suffer a loss $L>$ with probability $p>0$. The
% agent may save as much of his wealth s as he pleases, earning a
% per-period interest rate r. The agent may also insure any portion x of
% the uncertain loss the following period at a premium rate prem. That is,
% if the agent pays a premium prem*x this period, he receives an indemnity
% xL next period if he experiences a loss. The agent maximizes the present
% value of current and future utility of consumption, discounted at a
% subjective rate rho.
%
% States
%     w       stock of wealth
% Actions
%     s       savings
%     x       insurance coverage
% Parameters
%     alpha   relative risk aversion
%     rho     discount rate
%     r       interest rate on savings
%     ybar    income mean
%     L       magnitude of loss
%     p       probability of loss
%     load    premium load

function demdp21

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Numerical control parameters
n     = 500;                            % number of collocation nodes
maxit = 10;                             % maximum iterations in dpsolve
  
% Model parameters
smax  = 2.5;                            % maximum savings - 0 or 2.5-3.5
xmax  = 1;                              % maximum coverage - 0 or 1
alpha = 3.0;                           	% relative risk aversion
rho   = 0.10;                           % discount rate
r     = 0.08;                           % interest rate on savings
ybar  = 1.0;                            % income mean
L     = 0.5;                            % magnitude of insurable loss
p     = 0.1;                            % probability of insurable loss
load  = 0.2;                            % premium load
if r>rho, warning('Interest rate exceeds discount rate'), end

% Derived parameters
delta = 1/(1+rho);                      % discount factor
prem = (1+load)*delta*p*L;              % premium rate

% Model structure
model.func = @func;                     % model functions
model.params = {alpha ybar prem r smax xmax};% function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 2;                           % dimension of continuous action
model.e  = [0;L];                       % shocks
model.w  = [1-p;p];                     % shock probabilities

% Approximation structure
wmin  = ybar-L;                         % minimum state
wmax  = ybar+(1+r)*smax;                % maximum state
basis = fundefn('spli',n,wmin,wmax);    % basis functions
if wmin<0,error('Negative wealth in second period'), end

% Check model derivatives
dpcheck(model,1,[0 0])


%% SOLUTION

% Solve collocation equation
optset('dpsolve','maxit',maxit)
optset('dpsolve','algorithm','funcit')
[~,~,v,x] = dpsolve(model,basis);
optset('dpsolve','algorithm','newton')
[~,w,v,x,resid] = dpsolve(model,basis,v,x);


%% ANALYSIS

% Plot optimal savings and insurance policy
figure
hold on
plot(w,x(:,1)./w,'b','LineWidth',4)
plot(w,prem*x(:,2)./w,'r','LineWidth',4)
legend('Savings', 'Insurance','Location','NW')
title('Optimal Savings and Insurance Expenditures')
xlabel('Wealth')
ylabel('Expenditure as Proportion of Wealth')
xlim([wmin smax])

% Plot optimal insurance coverage policy
figure
plot(w,x(:,2),'r','LineWidth',4)
title('Optimal Insurance Coverage')
xlabel('Wealth')
ylabel('Coverage')
xlim([wmin smax])
ylim([0 1])

% Plot value function
figure
plot(w,v,'LineWidth',4)
title('Value Function')
xlabel('Wealth')
ylabel('Lifetime Utility')
xlim([wmin smax])

% Plot residual
figure
plot(w,resid,'LineWidth',4)
title('Residual')
xlabel('Wealth')
ylabel('Residual')
xlim([wmin smax])


%% SIMULATION

% Simulation parameters
nper = 10000;                           % number of periods simulated

% Initialize simulation
wtmp  = 1;                              % initial wealth
rng('default')                          
i = rand(nper,1)>1-p;

% Simulate model
wsim  = zeros(nper,1);
xsim  = zeros(nper,2);
for ip=1:nper
  wsim(ip) = wtmp;
  xsim(ip,:) = interp1(w,x,wtmp);
  wtmp = ybar + (1+r)*xsim(ip,1) - (1-xsim(ip,2))*i(ip)*L;
end
sse = [mean(xsim(:,1)./wsim); mean(prem*xsim(:,2)./wsim)];

% Print steady-state and ergodic moments
fprintf('\n\n')
fprintf('\n\n')
fprintf('ERGODIC MEAN EXPENDITURES AS PROPORTION OF WEALTH\n')
fprintf('   %-10s%8.3f\n','Savings     ',sse(1))
fprintf('   %-10s%8.3f\n','Insurance   ',sse(2))


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
        
function [out1,out2,out3] = func(flag,w,x,~,~,l,alpha,ybar,prem,r,smax,xmax)
n  = length(w);
switch flag
  case 'b'      % bounds
    out1 = zeros(n,2);
    out2 = [smax+zeros(n,1) xmax+zeros(n,1)];
    out3 = [];
  case 'f'      % reward
    out2 = zeros(n,2);
    out3 = zeros(n,2,2);
    c = w - x(:,1) - prem*x(:,2);
    dudc1 = c.^(-alpha);
    dudc2 = -alpha*c.^(-alpha-1);
    dcdx1 = -1;
    dcdx2 = -prem;
    out1        = (c.^(1-alpha))/(1-alpha);
    out2(:,1)   = dudc1*dcdx1;
    out2(:,2)   = dudc1*dcdx2;
    out3(:,1,1) = dudc2*dcdx1*dcdx1;
    out3(:,1,2) = dudc2*dcdx1*dcdx2;
    out3(:,2,1) = dudc2*dcdx2*dcdx1;
    out3(:,2,2) = dudc2*dcdx2*dcdx2;
  case 'g'      % transition
    out2 = zeros(n,1,2);
    out3 = zeros(n,1,2,2);
    out1 = ybar + (1+r)*x(:,1) - (1-x(:,2)).*l;
    out2(:,1,1) = (1+r)*ones(n,1);
    out2(:,1,2) = l;
end