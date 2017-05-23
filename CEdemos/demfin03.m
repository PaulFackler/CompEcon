% DEMFIN03 Heston's stochastic volatility option pricing model
function demfin03
close all
optset('finsolve','defaults');

disp(' ')
disp('DEMFIN03 Heston''s stochastic volatility option pricing model')

% Define parameters
r     = 0.05;     % risk free interest rate
delta = 0;        % rate of dividend payments
kappa = 1;        % mean-reversion parameter on volatility 
m     = 0.05;     % long-run mean volatility
sigma = 0.2;      % volatility of volatility
rho   = -0.5;     % correlation between price and volatility
K     = 1;        % strike price
put   = 1;        % 0=call, 1=put
T     = 1;        % time to maturity

% Create model variable
clear model
model.func='mffin03';
model.T=T;
model.american=0;
model.params={r,delta,kappa,m,sigma,rho,K,put};

% Define approximation space
n=[100 10];
smin=[log(0.01*K) 0.1*m];
smax=[log(5*K) 4*m];
fspace=fundefn('spli',n,smin,smax);
s=funnode(fspace);

% Call solution algorithm
N=fix(T*365+1);            % use one day time steps
N=1000;
c=finsolve(model,fspace,'implicit',s,N);

% Compute BS solution and produce plots
p=linspace(log(.5*K),log(1.5*K),251)';
nu=m;
S=exp(p);

V1=funeval(c(:,end),fspace,{p,nu});
V2=bs(sqrt(nu),S,K,r,delta,T,put);

close all
figure(1);
plot(S,V1,S,V2)
title('Option Values')
xlabel('S')
ylabel('Premium')
xlim([.75*K,1.5*K]);
legend('SV','Black-Scholes')

isigma=impvol(V1,S,K,r,delta,T,put,0);
figure(2);
plot(1./S,isigma,1./S,sqrt(nu)+zeros(size(S,1),1))
title('Implied Volatilities')
xlabel('K')
ylabel('Volatility')
xlim([0.75*K,1.5*K])
legend('SV','Black-Scholes')

nu=linspace(.1*m,4*m,7)';
nu=[.005 .05 .1 .125 .15 .175 .2]';
S=gridmake(p,nu);
V=reshape(funeval(c(:,end),fspace,S),251,7);
figure(3)
plot(exp(p),V)
title('Option Values for Alternative Values of \nu')
xlabel('S')
ylabel('Premium')
nn=length(nu);
legend([repmat('\nu = ',nn,1) reshape(sprintf('%5.3f',nu),5,nn)'],1)

prtfigs(mfilename,'Solution of the Stochastic Volatility Option Pricing Model',[1 2 3])

