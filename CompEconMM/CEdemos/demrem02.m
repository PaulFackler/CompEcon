%% DEMREM02 Commodity Storage Model
%
% Private expected profit maximizing storers enforce intertemporal price
% equilibriium in a market for a storable commodity.
%
% States
%     s         supply
% Response
%     x         ending stocks
% Parameters
%     gamma     inverse demand elasticity
%     xbar      storage capacity
%     k         unit storage cost
%     sigma     production volatility
%     delta     discount factor

function demrem02

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
gamma = 0.5;                                        % inverse demand elasticity
xbar  = 0.5;                                        % storage capacity
k     = 0.05;                                       % unit storage cost
sigma = 0.3;                                        % production  volatility
delta = 0.95;                                       % discount factor

% Continuous state shock distribution
m = 5;                                              % number of shocks
[y,w] = qnwlogn(m,-0.5*sigma^2,sigma^2);            % shocks and weights

% Model structure
model.func = @func;                                 % model functions
model.params = {delta,gamma,xbar,k};                % other parameters
model.ds = 1;                                       % dimension of state s
model.dx = 1;                                       % dimension of response x
model.e = y;                                        % shocks
model.w = w;                                        % shock probabilities

% Approximation structure
n      = 150;                                       % number of collocation nodes
smin   = 0.3;                                       % minimum state
smax   = 4.0;                                       % maximum state
[basis,Phi,s]  = fundefn('spli',n,smin,smax);       % basis functions


%% SOLUTION

% Solve rational expectations equilibrium directly by function iteration
tic
x = zeros(n,1);
c = zeros(n,1);
for it=1:100
  cold = c;
  f = -(s-x).^(-gamma)-k;
  d = -gamma*(s-x).^(-gamma-1);
  for j=1:m
    sn = x + y(j);
    f = f + w(j)*delta*funeval(c,basis,sn);
    d = d + w(j)*delta*funeval(c,basis,sn,1);
  end
  x = x + min(max(-f./d,-x),xbar-x);
  p = (s-x).^(-gamma);
  c = Phi\p;
  change = norm(c-cold,inf);
  fprintf ('%4i %10.1e\n',it,change)
  if change<1.e-10, break, end
end
toc

% Solve rational expectations equilibrium using remsolve
[c,s,x,f,resid] = remsolve(model,basis);
pe = (s-x).^(-gamma);
p0 = (s-0).^(-gamma);

% Critical supply levels
optset('broyden','output',0)
scrit1 = broyden(@(s) delta*w'*((   0+y-funeval(c,basis,   0+y)).^(-gamma))-(s-   0).^(-gamma)-k,1);
scrit2 = broyden(@(s) delta*w'*((xbar+y-funeval(c,basis,xbar+y)).^(-gamma))-(s-xbar).^(-gamma)-k,1);
pcrit1 = interp1(s,pe,scrit1);
pcrit2 = interp1(s,pe,scrit2);


%% ANALYSIS

% Plot equilibrium ending stocks function
figure
hold on
plot(s,x)
title('Equilibrium Ending Stocks')
xlabel('Supply')
ylabel('Ending Stocks')
xlim([0 smax])
ylim([0 xbar+0.1])

%  ... plot critical supply levels
plotvdash(scrit1,0)
plotvdash(scrit2,xbar)
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')

% Plot equilibrium price function
figure
hold on
plot(s,[pe p0])
legend('Total Demand','Consumption Demand')
title('Equilibrium Market Price')
xlabel('Supply')
ylabel('Price')
xlim([0 smax])

%  ... plot critical supply levels
plotvdash(scrit1,pcrit1)
plotvdash(scrit2,pcrit2)
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')

% Plot expected arbitrage profit function
figure
hold on
plot(s,f)
plothdash([],0)
title('Expected Arbitrage Profit')
xlabel('Supply')
ylabel('Profit')
xlim([0 smax])

%  ... plot critical supply levels
plotvdash(scrit1,0)
plotvdash(scrit2,0)
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')

% Plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
title('Response Residual')
xlabel('Supply')
ylabel('Residual')
xlim([0 smax])

%  ... plot critical supply levels
plotvdash(scrit1,0)
plotvdash(scrit2,0)
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')


%% SIMULATION

% Simulation parameters
nper = 20;                              % number of periods simulated
nrep = 10000;                           % number of replications

% Initialize simulation
sinit = ones(nrep,1);

% Simulate model directly, without remsimul
tic
ssim  = zeros(nrep,nper+1);
xsim  = zeros(nrep,nper+1);
rng('default')
ss = sinit;
for ip=1:nper+1
  xx = interp1(s,x,ss,'PCHIP');
  ssim(:,ip,:) = ss;
  xsim(:,ip,:) = xx;
  if ip<nper+1
    ss = xx + y(discrand(nrep,w),:);
  end
end
toc

% Simulate model using remsimul
tic
rng('default')
[ssim,xsim] = remsimul(model,basis,nper,sinit,s,x);
psim = (ssim-xsim).^-gamma;
toc

% Plot simulated and expected state path
figure
plot(0:nper,ssim(1:3,:),0:nper,mean(ssim),'k')
title('Simulated and Expected Supply')
xlabel('Period')
ylabel('Supply')

% Plot simulated and expected action path
figure
plot(0:nper,xsim(1:3,:),0:nper,mean(xsim),'k')
title('Simulated and Expected Ending Stocks')
xlabel('Period')
ylabel('Ending Stocks')

% Plot simulated and expected market price
figure
plot(0:nper,psim(1:3,:),0:nper,mean(psim),'k')
title('Simulated and Expected Market Price')
xlabel('Period')
ylabel('Market Price')

% Ergodic moments
savg = mean(ssim(:)); 
xavg = mean(xsim(:)); 
pavg = mean(psim(:)); 
sstd = std(ssim(:)); 
xstd = std(xsim(:)); 
pstd = std(psim(:)); 
fprintf('Ergodic Moments\n') 
fprintf('              Deterministic    Ergodic      Ergodic\n') 
fprintf('              Steady-State      Mean     Std Deviation\n') 
fprintf('Supply           %5.3f         %5.3f         %5.3f\n',[1 savg sstd])
fprintf('Ending Stocks    %5.3f         %5.3f         %5.3f\n',[0 xavg xstd])
fprintf('Market Price     %5.3f         %5.3f         %5.3f\n',[1 pavg pstd])

% Plot ergodic supply distribution
[qq,ss] = ksdensity(ssim(:),'support','positive','bandwidth',0.15);
figure
plot(ss,qq)
title('Ergodic Supply Distribution')
xlabel('Supply')
ylabel('Probability')
xlim([0 3])  

% Plot ergodic market price distribution
[qq,pp] = ksdensity(psim(:),'support','positive','bandwidth',0.15);
figure
plot(pp,qq)
title('Ergodic Market Price Distribution')
xlabel('Market Price')
ylabel('Probability')
xlim([0 3])  


%% SAVE FIGURES
printfigures(mfilename)


%% REMSOLVE FUNCTION FILE
%
%  User-supplied function that returns the bound, arbitrage, and state
%  transitions, and and their first derivatives with respect to the
%  response x, at an arbitrary number ns of states s and responses x
%  according to the format
%     [out1,out2,out3,out4] = func(flag,s,x,sn,xn,e,params)
%  where s is ns.ds states, x is ns.dx responses, sn is ns.ds states next
%  period, xn is ns.dx responses next period, and e is ns.de shocks.

function [out1,out2,out3,out4] = func(flag,s,x,sn,xn,y,delta,gamma,xbar,k)

ns = length(s);
switch flag
  case 'b'      % bounds
    out1 = zeros(ns,1);
    out2 = xbar*ones(ns,1);
  case 'f'      % arbitrage profit
    out1 = delta*(sn-xn).^(-gamma)-(s-x).^(-gamma)-k;
    out2 = -gamma*(s-x).^(-gamma-1);
    out3 = -gamma*delta*(sn-xn).^(-gamma-1);
    out4 =  gamma*delta*(sn-xn).^(-gamma-1);
  case 'g'      % transition
    out1 = x + y;
    out2 = ones(ns,1);
end