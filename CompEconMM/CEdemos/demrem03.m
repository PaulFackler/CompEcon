%% DEMREM03 Government Price Support Model
%
% A government buffer stock authority defends a minimum market price
% through the acquisition of excess supply.
%
% States
%     s         supply
% Responses
%     a         acreage planted
%     z         ending government stocks
% Parameters
%     zbar      maximum stocks
%     pbar      government support price
%     gamma     inverse consumption demand elasticity
%     beta      acreage supply elasticity
%     yvol      yield volatility
%     delta     discount factor

function demrem03

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
zbar  = 1.0;                                        % maximum stocks
pbar  = 0.9;                                        % government support price
gamma = 2.0;                                        % inverse consumption demand elasticity
beta  = 0.5;                                        % acreage supply elasticity
yvol  = 0.3;                                        % yield volatility
delta = 0.9;                                        % discount factor

% Continuous state shock distribution
m = 11;                                             % number of shocks
[y,w] = qnwlogn(m,-0.5*yvol^2,yvol^2);              % shocks and probabilities

% Model structure
model.func = @func;                                 % model functions
model.params = {delta zbar pbar gamma beta};        % function parameters
model.ds = 1;                                       % dimension of state s
model.dx = 2;                                       % dimension of response x
model.e = y;                                        % shocks
model.w = w;                                        % shock probabilities

% Approximation structure
n      = 501;                                       % number of collocation nodes
smin   = 0.2;                                       % minimum state
smax   = 5.0;                                       % maximum state
basis  = fundefn('spli',n,smin,smax);               % basis functions
s = funnode(basis);


%% SOLUTION

% Deterministic steady-state
pstar = broyden(@(pstar) (delta*pstar)^(1/beta)-pstar^(-1/gamma),1);
astar = (delta*pstar)^(1/beta);
if pstar>pbar, zstar=   0; end
if pstar<pbar, zstar=zbar; end
sstar = zstar+astar;

% Consumption demand at support price
dstar = pbar^(-1/gamma);

% Equilibrium government stocks and price
z = min(max(s-dstar,0),zbar);

% Solve rational expectations equilibrium directly by function iteration
maxit = 100;
a = ones(n,1);
for it=1:maxit
  Er = 0;
  Dr = 0;
  for j=1:m
    sn = z + a*y(j);
    zn = min(max(sn-dstar,0),zbar);
    pp = (sn-zn).^(-gamma);
    da = zeros(n,1);
    da(0>sn-dstar | sn-dstar>zbar) = y(j);
    dd = -gamma*((sn-zn).^(-gamma-1)).*da;
    Er = Er + w(j)*pp*y(j);
    Dr = Dr + w(j)*dd*y(j);
  end
  f = a - (delta*Er).^(1/beta);
  d = 1 - delta*(1/beta)*Dr.*(delta*Er).^((1/beta)-1);
  a = a - f./d;
  if norm(f)<1.e-10, break, end
end

% Solve rational expectations equilibrium using remsolve
x  = [0.5*ones(n,1) z];   
[c,s,x,f,resid] = remsolve(model,basis,x);
a  = x(:,1);
z  = x(:,2);
p = (s-z).^(-gamma);

% Critical supply levels
optset('broyden','output',0)
scrit1 = broyden(@(s) (s-funeval(c(:,2),basis,s)).^(-gamma)-pbar,1);
scrit2 = broyden(@(s) (s-zbar                   ).^(-gamma)-pbar,zbar+0.5);


%% ANALYSIS

% Plot equilibrium acreage planted
figure
hold on
plot(s,a)
title('Equilibrium Acreage Planted')
xlabel('Supply')
ylabel('Acreage')
xlim([0 smax])

%  ... plot critical supply levels
plotvdash(scrit1,max(a))
plotvdash(scrit2,min(a))
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')

% Plot equilibrium government ending stocks
figure
hold on
plot(s,z)
title('Equilibrium Ending Government Stocks')
xlabel('Supply')
ylabel('Stocks')
xlim([0 smax])  

%  ... plot critical supply levels
plotvdash(scrit1,0)
plotvdash(scrit2,zbar)
plottext(scrit1+0.02,[],'$s_1^*$')
plottext(scrit2+0.02,[],'$s_2^*$')

% Plot equilibrium price function
figure
hold on
plot(s,[p s.^(-gamma)])
legend('Total Demand','Consumption Demand')
plot([0 smax],[pbar pbar],'k:')
title('Equilibrium Market Price')
xlabel('Supply')
ylabel('Price')
xlim([0 smax])    
ylim([0 2])  

%  ... plot critical supply levels
plotvdash(scrit1,pbar)
plotvdash(scrit2,pbar)
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
nper = 30;                              % number of periods simulated
nrep = 10000;                           % number of replications

% Initialize simulation
sinit = ones(nrep,1);
rng('default')

% Preallocate arrays
ssim  = zeros(nrep,nper+1);
asim  = zeros(nrep,nper+1);
zsim  = zeros(nrep,nper+1);

% Simulate model
ss = sinit;
for ip=1:nper+1
  zz = min(max(ss-dstar,0),zbar);
  aa = interp1(s,a,ss,'PCHIP');
  ssim(:,ip,:) = ss;
  asim(:,ip,:) = aa;
  zsim(:,ip,:) = zz;
  if ip<nper+1
    ss = zz + aa.*y(discrand(nrep,w),:);
  end
end
psim = (ssim-zsim).^(-gamma);

% Plot simulated and expected state path
figure
plot(0:nper,ssim(1:2,:),0:nper,mean(ssim),'k')
title('Simulated and Expected Supply')
xlabel('Period')
ylabel('Supply')

% Plot simulated and expected response path 1
figure
plot(0:nper,asim(1:2,:),0:nper,mean(asim),'k')
title('Simulated and Expected Acreage Planted')
xlabel('Period')
ylabel('Acreage')

% Plot simulated and expected response path 2
figure
plot(0:nper,zsim(1:2,:),0:nper,mean(zsim),'k')
title('Simulated and Expected Ending Government Stocks')
xlabel('Period')
ylabel('Ending Stocks')

% Plot simulated and expected market price
figure
plot(0:nper,psim(1:2,:),0:nper,mean(psim),'k')
title('Simulated and Expected Market Price')
xlabel('Period')
ylabel('Market Price')

% Ergodic moments
savg = mean(ssim(:));
aavg = mean(asim(:));
zavg = mean(zsim(:));
pavg = mean(psim(:));
sstd = std(ssim(:));
astd = std(asim(:));
zstd = std(zsim(:));
pstd = std(psim(:));
fprintf('Ergodic Moments\n') 
fprintf('              Deterministic    Ergodic      Ergodic\n') 
fprintf('              Steady-State      Mean     Std Deviation\n') 
fprintf('Supply           %5.3f         %5.3f         %5.3f\n',[sstar savg sstd])
fprintf('Acreage Planted  %5.3f         %5.3f         %5.3f\n',[astar aavg astd])
fprintf('Ending Stocks    %5.3f         %5.3f         %5.3f\n',[zstar zavg zstd])
fprintf('Market Price     %5.3f         %5.3f         %5.3f\n',[pstar pavg pstd])

% Plot ergodic state distribution
[qq,ss] = ksdensity(ssim(:),'support','positive','bandwidth',0.14);
figure
plot(ss,qq)
title('Ergodic Supply Distribution')
xlabel('Supply')
ylabel('Probability')
xlim([0 4])  


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

function [out1,out2,out3,out4] = func(flag,s,x,sn,xn,e,delta,zbar,pbar,gamma,beta)
a  = x(:,1);
z  = x(:,2);
ns = length(s);
ds = 1;
dx = 2;
switch flag
  case 'b'      % bounds
    out1 = zeros(ns,2);
    out2 = [inf+zeros(ns,1) zbar+zeros(ns,1)];
  case 'f'      % arbitrage profit
    zn = xn(:,2);
    out1 = zeros(ns,dx);
    out2 = zeros(ns,dx,dx);
    out3 = zeros(ns,dx,ds);
    out4 = zeros(ns,dx,dx);
    p     = (s -z ).^(-gamma);
    pn    = (sn-zn).^(-gamma);
    pder  = -gamma*(s -z ).^(-gamma-1);
    pnder = -gamma*(sn-zn).^(-gamma-1);
    out1(:,1)   = delta*pn.*(sn-z)./a-a.^beta;
    out1(:,2)   = pbar-p;
    out2(:,1,1) = -delta*pn.*(sn-z)./(a.^2)-beta*a.^(beta-1);
    out2(:,1,2) = -delta*pn./a;
    out2(:,2,2) = pder;
    out3(:,1,1) = delta*(pn+pnder.*(sn-z))./a;
    out4(:,1,2) = -delta*pnder.*(sn-z)./a;
  case 'g'      % transition
    out2 = zeros(ns,ds,dx);
    out1 = z+a.*e;
    out2(:,1,1) = e;
    out2(:,1,2) = ones(ns,1);
end