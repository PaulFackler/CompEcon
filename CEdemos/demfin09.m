% DEMFIN09 Financial asset calibration demonstration
close all
optset('finsolve','defaults');

disp(' ')
disp('DEMFIN09 Financial asset calibration demonstration')
close all

% Define parameters
kappa = 0.0363;
alpha = 0.0692; 
sigma = 0.0272;
 
% Create model variable
clear model
model.func='mffin01';
model.T=30;
model.american=0;
model.params={kappa,alpha,sigma};

% Define approximation space
n=20;
fspace=fundefn('cheb',n,0,2);
s=funnode(fspace);

% Call solution algorithm
optset('finsolve','keepall',1);
c=finsolve(model,fspace,'lines',s,120);

% Calibrate to the data
y=[4.44	 4.49 	4.51	4.63	4.63	4.62	4.82	4.77	5.23];
tau=[.25 .5 1 2 3 5 7 10 30];
V=exp(-y/100.*tau);
m=length(V);
t=(0:0.25:30);
tind=[2 3 5 9 13 21 29 41 121];
% t(tind)=tau and columns of c correspond to t
s=findstate(c(:,tind),fspace,alpha,V);

% Create plots
Vhat=funeval(c,fspace,s);
warning off
yhat=-100*log(Vhat)./t; yhat(1)=100*s;
warning on

figure(1)
plot(tau,y,'*',t,yhat)
title('Actual and Fitted Bond Yields')
xlabel('Time to Maturity')
ylabel('Yield')

disp('Model short interest rate and 3-month rate')
disp([s*100 y(1)])

prtfigs(mfilename,'Actual and Fitted Bond Yields',[1])