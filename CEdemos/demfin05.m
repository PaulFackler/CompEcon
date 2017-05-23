% DEMFIN05 Barrier option pricing problem
function demfin05
close all
optset('finsolve','defaults');

disp(' ')
disp('DEMFIN05 Barrier option pricing problem')

% Define parameters
r      = 0.05;               % risk free interest rate
delta  = 0;                  % dividend rate 
sigma  = 0.2;                % volatility  
K      = 1;                  % exercise price 
put    = 0;                  % put indicator
T      = 1;                  % time-to-maturity
Sb     = 0.8;                % barrier 

% Create model variable
clear model
model.func='mffin04';
model.T=T;
model.american=0;
model.params={r,delta,sigma,K,put};
model.barrier=[0 Sb 1;0 inf 1];

% Define approximation space
n=75; 
fspace=fundefn('spli',n,Sb,2*K);    
S=funnode(fspace);

% Call solution algorithm
[c,Vhat]=finsolve(model,fspace);

% Create plots
S=sort([linspace(.5*K,2*K,501)';Sb-eps;Sb]);
Vhat=funeval(c,fspace,S);
Vhat(S<Sb)=0;
Vbs=bs(sigma,S,K,r,delta,T,put);
V=Vbs-(S./Sb).^(1-2*r/sigma^2).*bs(sigma,Sb^2./S,K,r,delta,T,put);
V(S<Sb)=0;

figure(1)
plot(S,Vhat,S,Vbs);
title('Down-and-out Call Option Premium')
xlabel('S')
ylabel('Pemium')
legend({'Barrier Option','Vanilla Option'},2);
xlim([.5 1.5]*K)
ylim([0 0.5])

figure(2)
plot(S,[V-Vhat])
title('Approximation Error')
xlabel('S')
ylabel('Error')
xlim([Sb,2*K])

disp('Sb Vb(Sb) Vbs(Sb)')
disp([S(S==Sb) Vhat(S==Sb) Vbs(S==Sb)])

n=75; 
lambda=.5;
Delta=(2*K-Sb)/(n-2);
fspace=fundefn('lin',n,Sb-lambda*Delta,2*K);    
% Call solution algorithm
s=funnode(fspace);
s(abs(s-Sb)<1e-13)=Sb;
[c,Vhat]=finsolve(model,fspace,[],s);

Vhat=funeval(c,fspace,S);
figure(3)
plot(S,[V-Vhat])
title('Approximation Error with Barrier Not a Node')
xlabel('S')
ylabel('Error')
xlim([Sb,2*K])

prtfigs(mfilename,'Solution of the Barrier Option Pricing Model',[1 2 3])


