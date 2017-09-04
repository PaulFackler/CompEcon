%% DEMDP22 Savings and Insurance - Sensitivity Analysis
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

function demdp22

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
smax  = 3.5;                            % maximum savings - 0 or 2.5-3.5
xmax  = 1;                              % maximum coverage - 0 or 1
alpha = 3.0;                           	% relative risk aversion
rho   = 0.10;                           % discount rate
r     = 0.08;                           % interest rate on savings
ybar  = 1.0;                            % income mean
L     = 0.5;                            % magnitude of insurable loss
p     = 0.2;                            % probability of insurable loss
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
n     = 200;                            % number of collocation nodes
wmin  = ybar-L;                         % minimum state
wmax  = ybar+(1+r)*smax;                % maximum state
basis = fundefn('spli',n,wmin,wmax);    % basis functions
if wmin<0,error('Negative wealth in second period'), end

% Check model derivatives
dpcheck(model,1,[0 0])

% Numerical control parameters
nper  = 1000;                           % number of periods simulated
xlm   = [0 2.5];                        % x-axis plot limits
ylm   = [0 0.7];                        % y-axis plot limits


%% SOLUTION

% Both savings and insurance available
[~,wA,vA,xA,residA,sseA] = solve(ybar,L,p,model,basis,nper,r,prem);

% Only insurance available, no savings
model.params = {alpha ybar prem r 0 xmax};
[~,wI,vI,xI,residI,sseI] = solve(ybar,L,p,model,basis,nper,r,prem);

% Only savings available, no insurance
model.params = {alpha ybar prem r smax 0};
[~,wS,vS,xS,residS,sseS] = solve(ybar,L,p,model,basis,nper,r,prem);

% Neither savings nor insurance available
model.params = {alpha ybar prem r 0 0};
[~,w0,v0,~,resid0,~] = solve(ybar,L,p,model,basis,nper,r,prem);


%% OUTPUT

% Print ergodic means of savings and insurance expenditures
fprintf('\n\n')
fprintf('\n\n')
fprintf('ERGODIC MEAN EXPENDITURES AS PROPORTION OF WEALTH\n')
fprintf('               Savings and  Insurance    Savings\n')
fprintf('                Insurance     Only        Only\n')
fprintf('%-10s%8.3f%11.3f%12.3f\n','Savings        ',sseA(1),sseI(1),sseS(1))
fprintf('%-10s%8.3f%11.3f%12.3f\n','Insurance      ',sseA(2),sseI(2),sseS(2))

% Plot optimal savings and insurance policy
figure
hold on
plot(wA,xA(:,1)./wA,'b','LineWidth',4)
plot(wA,prem*xA(:,2)./wA,'r','LineWidth',4)
plot(wS,xS(:,1)./wS,'b:','LineWidth',4)
plot(wI,prem*xI(:,2)./wI,'r:','LineWidth',4)
legend('Savings', 'Insurance', ...
       'Savings w/o Insurance','Insurance w/o Savings','Location','NW')
xlim(xlm)
ylim(ylm)
title('Optimal Savings and Insurance Expenditures')
xlabel('Wealth')
ylabel('Expenditure as Proportion of Wealth')

% Plot optimal insurance coverage policy
figure
hold on
plot(wA,xA(:,2),'r','LineWidth',4)
plot(wI,xI(:,2),'r:','LineWidth',4)
legend('With Savings', 'Without Savings','Location','NW')
xlim(xlm)
ylim([0 1])
title('Optimal Insurance Coverage')
xlabel('Wealth')
ylabel('Coverage')

% Plot value function
figure
plot(wA,vA,wS,vS,wI,vI,w0,v0,'LineWidth',4)
legend('Savings and Insurance','Savings Only','Insurance Only','Neither','Location','SE')
xlim(xlm)
title('Value Function')
xlabel('Wealth')
ylabel('Lifetime Utility')

% Plot residuals
figure
plot(wA,residA,wS,residS,wI,residI,w0,resid0,'LineWidth',4)
legend('Savings and Insurance','Savings Only','Insurance Only','Neither','Location','SE')
xlim(xlm)
title('Residual')
xlabel('Wealth')
ylabel('Residual')


%% PARAMETRIC ANALYSIS

% Number of parameters simulated
nsim  = 15;

% Ergodic mean savings and insurance expenditures vs. premium load
loadmin = -0.5;
loadmax =  0.5;
loadsim = nodeunif(nsim,loadmin,loadmax);
covsim = zeros(1,nsim);
ssesim = zeros(2,nsim);
for isim=1:nsim
  prem = (1+loadsim(isim))*delta*p*L;
  model.params = {alpha ybar prem r smax xmax};
  [~,~,~,~,~,ssesim(:,isim),covsim(isim)] = solve(ybar,L,p,model,basis,nper,r,prem);
end

% Plot ergodic mean savings and insurance expenditures vs. premium load
figure
plot(loadsim,ssesim,'LineWidth',4)
legend('Savings','Insurance','Location','E')
title('Ergodic Mean Savings and Insurance Expenditures vs. Premium Load')
xlabel('Premium Load')
ylabel('Expenditure as Proportion of Wealth')

% Plot ergodic mean coverage vs. premium load
figure
plot(loadsim,covsim,'LineWidth',4)
title('Ergodic Mean Coverage vs. Premium Load')
xlabel('Premium Load')
ylabel('Proportion of Loss Insured')

% Ergodic mean savings and insurance expenditures vs. interest rate
rmin = 0.00;
rmax = 0.08;
rsim = nodeunif(nsim,rmin,rmax);
covsim = zeros(1,nsim);
ssesim = zeros(2,nsim);
for isim=1:nsim
  model.params = {alpha ybar prem rsim(isim) smax xmax};
  [~,~,~,~,~,ssesim(:,isim),covsim(isim)] = solve(ybar,L,p,model,basis,nper,rsim(isim),prem);
end

% Plot ergodic mean savings and insurance expenditures vs. interest rate
figure
plot(rsim,ssesim,'LineWidth',4)
legend('Savings','Insurance','Location','E')
title('Ergodic Mean Savings and Insurance Expenditures vs. Interest Rate')
xlabel('Interest Rate')
ylabel('Expenditure as Proportion of Wealth')

% Plot ergodic mean coverage vs. interest rate
figure
plot(rsim,covsim,'LineWidth',4)
title('Ergodic Mean Coverage vs. Interest Rate')
xlabel('Interest Rate')
ylabel('Proportion of Loss Insured')


%% SAVE FIGURES
printfigures(mfilename)


%% ANCILLARY FUNCTION
%
% Solves and simulates model

function [c,w,v,x,resid,sse,cov] = solve(ybar,L,p,model,basis,nper,r,prem)

% Solve collocation equation
optset('dpsolve','output',0)
optset('dpsolve','maxit',15)
optset('dpsolve','algorithm','funcit')
[~,~,v,x] = dpsolve(model,basis);
optset('dpsolve','algorithm','newton')
[c,w,v,x,resid] = dpsolve(model,basis,v,x);

% Simulate model
if nper>0
  rng('default')
  wsim  = zeros(nper,1);
  xsim  = zeros(nper,2);
  i = rand(nper,1)>1-p;
  wtmp  = 1;
  for ip=1:nper
    wsim(ip) = wtmp;
    xsim(ip,:) = interp1(w,x,wtmp);
    wtmp = ybar + (1+r)*xsim(ip,1) - (1-xsim(ip,2))*i(ip)*L;
  end
  sse = [mean(xsim(:,1)./wsim); mean(prem*xsim(:,2)./wsim)];
  cov = mean(xsim(:,2));
end


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