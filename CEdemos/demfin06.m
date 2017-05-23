% DEMFIN06 Compound option pricing example - option on a bond
function demfin07
close all
optset('finsolve','defaults');

disp(' ')
disp('DEMFIN06 Compound option pricing example - option on a bond')

% Define parameters 
kappa = 0.1;       % speed of mean reversion
alpha = 0.05;      % long-run mean interest rate
sigma = 0.1;       % interest rate volatility 
TB    = 30;        % bond maturity
K     = 0.2;       % exercise price
put   = 0;         % put indicator
TO    = 1;         % option maturity

% Create model variable for the bond
clear model
modelB.func='mffin06';
modelB.T=TB;
modelB.american=0;
modelB.params={kappa,alpha,sigma};

% Define approximation space for the bond
n=20;
fspaceB=fundefn('cheb',n,0,2);

% Call the solver
cB=finsolve(modelB,fspaceB);

% Create model variable the option
modelO.func='mffin06';
modelO.T=TO;
modelO.american=0;
modelO.params={kappa,alpha,sigma,K,put,cB,fspaceB};

% Define approximation space for the option
n=80;
fspaceO=fundefn('spli',n,0,2);

% Call the solver
cO=finsolve(modelO,fspaceO);

% Create plots
x=linspace(0,.25,201)';
Vhat=funeval(cO,fspaceO,x);

figure(1)
plot(x,Vhat); 
title(['As a Function of Short Rate'])
xlabel('r')
ylabel('Option Premium')
ylim([0 0.25])

BondVal=funeval(cB,fspaceB,x);
figure(2)
plot(BondVal,Vhat); 
title(['As a Function of Bond Value'])
xlabel('Bond Value')
ylabel('Option Premium')
ylim([0 0.25])
xlim([0 .45])

prtfigs(mfilename,'Solution of the Compound Option Pricing Model',[1 2])