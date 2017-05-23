% DEMFIN07 Asian option pricing demonstration
close all
optset('finsolve','defaults');

disp(' ')
disp('DEMFIN07a Asian option pricing demonstration')

% Define parameters
r     = 0.1;
delta = 0;
sigma = 0.4;
L     = 1/2;
put   = 1;
tau   = 1/4;

% Create model variable
clear model
model.func='mffin07a';
model.T=tau;
model.american=1;
model.params={r,delta,sigma,L,put};

% Define approximation space
n=21; 
fspace=fundefn('spli',[101 101],[0.001 0.001],[5 2.5]);

% Call solution algorithm
c=finsolve(model,fspace,'implicit',[],101);

% Transform solution to natural state space (y to S)
S=linspace(0.001,2,101)';
A=[0.5 0.75 1 1.25 1.5]*((L-tau)/L);
m=length(A);
Vhat=zeros(length(S),m);
for i=1:m
  Vhat(:,i)=funeval(c,fspace,{S,A(i)});
%  if ~put, Vhat(y>fspace.b,i)=0; end
end

% Create plots
figure(1)
plot(S,Vhat);
title('Call Option Premium')
xlabel('S')
ylabel(['V(S,A,' num2str(tau) ')']);
nn=length(A);
legend([repmat('M = ',nn,1) reshape(sprintf('%4.2f',A*(L/(L-tau))),4,nn)'],2)
set(gca,'ylim',max(0,get(gca,'ylim')));

