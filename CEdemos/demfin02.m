% DEMFIN02 Black-Scholes option pricing model using function approximation
function demfin02
close all
optset('finsolve','defaults');

disp(' ')
disp('DEMFIN02 Black-Scholes option pricing model using function approximation')

% Define parameters
r      = 0.05;               % risk free interest rate
deltaS = 0;                  % dividend rate 
sigma  = 0.2;                % volatility  
K      = 1;                  % exercise price 
put    = 0;                  % put indicator
T      = 1;                  % time-to-maturity

% Create model variable
clear model
model.func='mffin02';
model.T=T;
model.params={r,deltaS,sigma,K,put};
model.american=0;

% Define approximation space
n=51; 
fspace=fundefn('lin',n,0,2*K);    
s=funnode(fspace);

% Call solution algorithm
c=finsolve(model,fspace,'lines',s);

% Compute exact solution and plot solution and errors
S=linspace(0,2*K,501)';
premium=bs(sigma,S,K,r,deltaS,T,put);

figure(1)
plot(S,funeval(c,fspace,S));
title('Call Option Premium')
xlabel('S')
ylabel('Premium')

figure(2)
plot(S,[premium-funeval(c,fspace,S)])
title('Approximation Errors')
xlabel('S')
ylabel('Error')

% Compute performance comparisons

stype={'lines' 'explicit' 'implicit' 'CN' 'stiff'};
N=[1 75 75 75 75];
m=length(stype);
C=zeros(n,m);
tl=zeros(1,m);

fspace=fundefn('lin',n,0,2*K);    
s=funnode(fspace);

for i=1:m
  tic;
  c=finsolve(model,fspace,stype{i},s,N(i));
  tl(i)=toc;
  C(:,i)=c(:,end);
end
VL=funeval(C,fspace,S);

ts=zeros(1,m);
N=[1 250 75 75 75];

fspace=fundefn('spli',n,0,2*K);    
s=funnode(fspace);
for i=1:m
  tic;
  c=finsolve(model,fspace,stype{i},s,N(i));
  ts(i)=toc;
  C(:,i)=c(:,end);
end
VS=funeval(C,fspace,S);

disp('Maximum errors on [0,2K] for lines, explicit, implicit, CN and stiff methods')
disp('Piecewise Linear')
show(max(abs(premium(:,ones(5,1))-VL)),6)
disp('Cubic Spline')
show(max(abs(premium(:,ones(5,1))-VS)),6)

disp('Timing for lines, explicit, implicit, CN and stiff methods')
disp('Piecewise Linear')
show(tl)
disp('Cubic Spline')
show(ts)

prtfigs(mfilename,'Black-Scholes Option Pricing Model: Approximation Errors',[2])
