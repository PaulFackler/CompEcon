% DEMFIN07 Asian option pricing demonstration
%close all
%optset('finsolve','defaults');

disp(' ')
disp('DEMFIN07 Asian option pricing demonstration')

% Define parameters
r     = 0.1;
delta = 0;
sigma = 0.15;
L     = 1/2;
put   = 1;
tau   = 1/4;


tau=1/4;


% Create model variable
clear model
model.func='mffin07';
model.T=tau;
model.american=0;
model.params={r,delta,sigma,L,put};

% Define approximation space
n=201; 
fspace=fundefn('lin',n,0,1);

% Call solution algorithm
c=finsolve(model,fspace);

M=[.5 .75 1 1.25 1.5];
funeval(c,fspace,(L-tau)*M')'
return
% Transform solution to natural state space (y to S)
S=linspace(0.001,2,101)';
M=[.5 .75 1 1.25 1.5];
m=length(M);
Vhat=zeros(length(S),m);
for i=1:m
  y=((L-tau)*M(i))./S;
  Vhat(:,i)=S.*funeval(c,fspace,y);
  if ~put, Vhat(y>fspace.b,i)=0; end
end

% Create plots
figure(1)
plot(S,Vhat);
title('Call Option Premium')
xlabel('S')
ylabel(['V(S,M,' num2str(tau) ')']);
nn=length(M);
legend([repmat('M = ',nn,1) reshape(sprintf('%4.2f',M),4,nn)'],2)
set(gca,'ylim',max(0,get(gca,'ylim')));

y=linspace(fspace.a,fspace.b,101)';
figure(2)
plot(y,funeval(c,fspace,y));
title('Approximation Function')
xlabel('y')
ylabel('v(y,\tau)');
ylim([-0.1 0.6])

%prtfigs(mfilename,'Solution of the Asian Option Pricing Model',[2 1])
