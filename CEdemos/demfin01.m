% DEMFIN01 Cox-Ingersoll-Ross bond pricing example
%function demfin01
%close all

disp(' ')
disp('DEMFIN01 Cox-Ingersoll-Ross bond pricing example')

% Define parameters
T=30;
kappa=.1;
alpha=.05;
sigma=0.1;

% Evaluate for values of the short rate between 0 and 2
%   (2 is as close to infinity as one needs to get for interest rates
%      unless the economy is in hyperinflation.)
n=20;
fspace=fundefn('cheb',n,0,2);

if 1
  % Use finsolve
  clear model
  model.func='mffin01';
  model.T=T;
  model.american=0;
  model.params={kappa,alpha,sigma};
  s=funnode(fspace);
  c=finsolve(model,fspace,'lines',s,1);
else
  % use the specialized function cirbond
  c=cirbond(fspace,T,kappa,alpha,sigma);
end

% Compute exact solution and plot solution and errors
sigma2=sigma*sigma;
gam=sqrt(kappa.^2+2*sigma2);
egam=exp(gam*T)-1;
denom=(gam+kappa).*egam+2*gam;
B=2*egam./denom;
A=(2*gam*exp((kappa+gam)*T/2)./denom).^(2*kappa*alpha./sigma2);

x=linspace(0,0.25,101)';
V=A*exp(-B*x);
Vhat=funeval(c(:,end),fspace,x);

figure(1)
plot(x,Vhat); 
title('Bond Price')
xlabel('r')
ylabel('V(r)')

figure(2)
plot(x,V-Vhat)
title('Approximation Errors')
xlabel('r')
ylabel('Error')

fprintf('Maximum error on [0,0.25]: %10.4e\n',max(abs(V-Vhat)))

prtfigs(mfilename,['Solution of the CIR ' num2str(T) ' Year Bond Pricing Model'],[1 2])