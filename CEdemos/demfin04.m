% DEMFIN04 American put option pricing problem
% using function approximation
function demfin04
close all
optset('finsolve','defaults');

disp(' ')
disp('DEMFIN04 American put option pricing problem')

% Define parameters
r      = 0.05;               % risk free interest rate
deltaS = 0;                  % dividend rate 
sigma  = 0.2;                % volatility  
K      = 1;                  % exercise price 
put    = 1;                  % put indicator
T      = 1;                  % time-to-maturity

% Create model variable
clear model
model.func='mffin04';
model.T=T;
model.american=1;
model.params={r,deltaS,sigma,K,put};

% Define approximation space
n=100; 
fspace=fundefn('lin',n,0,2*K);    
s=funnode(fspace);

% Call solution algorithm
N=75;
optset('finsolve','keepall',1);
c=finsolve(model,fspace,[],s,N);

% Compute Barone-Adesi/Whaley solution and plot solution and errors
S=linspace(fspace.a,fspace.b,301)';
V=funeval(c(:,end),fspace,S);
[Vbaw,bs,sstarB]=baw(sigma,S,K,r,deltaS,T,put);

figure(1)
plot(S,V,'k-',S,bs,'r-',sstarB,funeval(c(:,end),fspace,sstarB),'k*');
title('Option Premium')
xlabel('S')
ylabel('Premium')
legend('American','European')

figure(2)
plot(S,[Vbaw-V],'k')
title('Approximation Error')
xlabel('S')
ylabel('Error')

% Compute and plot the optimal exercise boundary
V0=feval(model.func,'V0',s,model.params{:});
temp=funeval(c,fspace,s)==V0(:,ones(1,N+1)) & s(:,ones(1,N+1))<=K;
% use upper bound on S* at time values where the upper bound changes
sstar=s(sum(temp)'+1);
sstarl=s(sum(temp)');
tau=linspace(0,T,N+1)';
ind=[1;find(diff(sstar)~=0)+1;N+1];
sstari=sstar(ind);
taui=tau(ind);
% end point adjustments
sstari(1)=K;
send=sstari(end-1)+(sstari(end-1)-sstari(end-2))/(taui(end-1)-taui(end-2))*(taui(end)-taui(end-1));
sstari(end)=(sstari(end)+send)/2;

figure(3)
plot(taui,sstari,'k')
hold on; stairs(tau,sstar,'r--'); stairs(tau,sstarl,'r--'); hold off
title('Early Exercise Boundary')
xlabel('\tau')
ylabel('S^*')
ylim([0.8 1.05])

prtfigs(mfilename,'Solution of the American Put Option Pricing Model',[1 2 3])