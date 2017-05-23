% DEMBVP02 Commodity storage example
% X=[p;s]
%   p: price
%   s: stocks
%   p'(t)=r*p(t)+C
%   s'(t)=-D(p(t)) = -p(t)^(-eta)
% s.t.
%   s(0)=S0
%   s(1)=0
function dembvp02
 disp(' ')
 disp('DEMBVP02 Commodity storage example')

% Set parameter values
 r   = 0.10;
 C   = 0.5;
 eta = 2;
 S0  = 1;
 T   = 1;

% Define model structure
 clear model
 model.func='pbvp02';
 model.tb=[0;T];
 model.params={r,C,eta,S0};

% Define approximation space
 n=6;
 fspace=fundefn('cheb',n-1,0,T);
 tnodes=funnode(fspace);
 fspace=fundefn('cheb',n,0,T);

 N=501;
 t=linspace(0,T,N)';
% Initial conditions: constant price and consumption
 cc=funfitxy(fspace,t,[(S0/T).^(-1/eta)+zeros(N,1) 1-(S0/T)*t]);

% Call solver
 optset('broyden','defaults')
 [c,x,r]=bvpsolve(model,fspace,tnodes,cc,t);

% Create plots
 close all
 figure(1);
 plot(t,x)
 title('Equilibrium Price and Stock Level')
 xlabel('t')
 text(0.68,1.1,'price')
 text(0.5,0.4,'stocks')

 figure(2)
 plot(t,r(:,1)*1e6,'-',t,r(:,2),'--')
 title('Residual Functions')
 legend('P'' x 10^6','S''')
 xlabel('t')

% Display residual summary information
 disp('Maximum absolute residuals for price and stock level')
 disp([max(abs(r(:,1))) max(abs(r(:,2)))])
