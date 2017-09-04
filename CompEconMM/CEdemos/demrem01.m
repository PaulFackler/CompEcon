%% DEMREM01 Asset Pricing Model
%
% Lucas-Prescott rational expectations asset pricing model.
%
% States
%     d       asset dividend
% Response
%     p       asset price
% Parameters
%     beta    coefficient of risk aversion
%     dbar    long-run mean dividend
%     gamma   dividend autoregression coefficient
%     sigma   dividend shock standard deviation
%     delta   discount factor

function demrem01

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
beta  = 0.5;                                        % coefficient of risk aversion
dbar  = 1.0;                                        % long-run mean dividend
gamma = 0.5;                                        % dividend autoregression coefficient
sigma = 0.1;                                        % dividend shock standard deviation
delta = 0.9;                                        % discount factor

% Continuous state shock distribution
m = 5;                                              % number of shocks
[e,w] = qnwnorm(m,0,sigma^2);                       % shocks and weights

% Model structure
model.func = @func;                                 % model functions
model.params = {delta beta dbar gamma};             % function parameters
model.ds = 1;                                       % dimension of state s
model.dx = 1;                                       % dimension of response x
model.e  = e;                                       % shocks
model.w  = w;                                       % shock probabilities

% Approximation structure
n      = 25;                                        % number of collocation nodes
dmin   = 0.1;                                       % minimum state
dmax   = 1.9;                                       % maximum state
[basis,Phi,dnode]  = fundefn('cheb',n,dmin,dmax);   % basis functions


%% SOLUTION

% Deterministic steady-state
pstar = delta*dbar/(1-delta);                       % asset price

% Solve rational expectations equilibrium directly
tic
LHS = diag(dnode.^(-beta))*Phi;
RHS = 0;
for k=1:m
  dnext = dbar + gamma*(dnode-dbar) + e(k);
  LHS   = LHS - delta*w(k)*diag(dnext.^(-beta))*funbase(basis,dnext);
  RHS   = RHS + delta*w(k)*dnext.^(1-beta);
end
c = LHS\RHS;
p = funeval(c,basis,dnode);

% Residual
dd = nodeunif(10*n,dmin,dmax);
Ef = 0;
for k=1:m
  dnext = dbar + gamma*(dd-dbar) + e(k);
  f     = diag(dnext.^(-beta))*(funeval(c,basis,dnext)+dnext);
  Ef    = Ef + delta*w(k)*f;
end
resid = dd.^(-beta).*funeval(c,basis,dd)-Ef;
toc

% Solve rational expectations equilibrium using remsolve
[c,d,p,f,resid] = remsolve(model,basis,p);


%% ANALYSIS

% Plot equilibrium asset price
figure
plot(d,p)
title('Equilibrium Asset Price')
xlabel('Dividend')
ylabel('Price')
xlim([0.1 1.9])   

% Plot expected arbitrage profit function
figure
hold on
plot(d,f)
plothdash([],0)
title('Expected Arbitrage Benefit')
xlabel('Dividend')
ylabel('Profit')
xlim([0.1 1.9])   

% Plot residual
figure
hold on
plot(d,resid)
plothdash([],0)
title('Response Residual')
xlabel('Dividend')
ylabel('Residual')
xlim([0.1 1.9])   


%% SIMULATION

% Simulation parameters
nper = 50;                              % number of periods simulated
nrep = 10000;                           % number of replications

% Initialize simulation
dinit = ones(nrep,1);
rng('default')

% Preallocate arrays
dsim = zeros(nrep,nper);
psim = zeros(nrep,nper);

% Simulate model directly, without remsimul
dd = dinit;
for ip=1:nper+1
  pp = interp1(d,p,dd,'PCHIP');
  dsim(:,ip,:) = dd;
  psim(:,ip,:) = pp;
  if ip<nper+1
    dd = dbar+gamma*(dd-dbar)+e(discrand(nrep,w),:);
  end
end

% Simulate model using remsimul
rng('default')
[dsim,psim] = remsimul(model,basis,nper,dinit,d,p);

% Plot simulated and expected state path
figure
plot(0:nper,dsim(1:3,:),0:nper,mean(dsim),'k')
title('Simulated and Expected Asset Dividend')
xlabel('Period')
ylabel('Dividend')

% Plot simulated and expected response path
figure
plot(0:nper,psim(1:3,:),0:nper,mean(psim),'k')
title('Simulated and Expected Asset Price')
xlabel('Period')
ylabel('Price')

% Ergodic moments
davg = mean(dsim(:)); 
pavg = mean(psim(:)); 
dstd = std(dsim(:)); 
pstd = std(psim(:)); 
fprintf('Ergodic Moments\n') 
fprintf('              Deterministic    Ergodic      Ergodic\n') 
fprintf('              Steady-State      Mean     Std Deviation\n') 
fprintf('Asset Dividend   %5.3f         %5.3f         %5.3f\n',[dbar  davg dstd])
fprintf('Asset Price      %5.3f         %5.3f         %5.3f\n',[pstar pavg pstd])

% Plot ergodic dividend distribution
[qq,dd] = ksdensity(dsim(:),'support','positive','bandwidth',0.04);
figure
plot(dd,qq)
title('Ergodic Dividend Distribution')
xlabel('Dividend')
ylabel('Probability')
xlim([0 2])   

% Plot ergodic asset price distribution
[qq,pp] = ksdensity(psim(:),'support','positive','bandwidth',0.04);
figure
plot(pp,qq)
title('Ergodic Asset Price Distribution')
xlabel('Price')
ylabel('Probability')
xlim([5 15])  


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

function [out1,out2,out3,out4] = func(flag,d,p,dn,pn,e,delta,beta,dbar,gamma)

ns = length(d);
switch flag
  case 'b'      % bounds
    out1 = zeros(ns,1)-inf;
    out2 = zeros(ns,1)+inf;
  case 'f'      % reward
    u = d.^(-beta);
    un = dn.^(-beta);
    out1 = p.*u-delta*(pn+dn).*un;
    out2 = u;
    out3 = beta*delta*(pn+dn).*un./dn-delta*un;
    out4 = -delta*un;
  case 'g'      % transition
    out1 = dbar+gamma*(d-dbar)+e;
    out2 = zeros(ns,1);
end