% DEMAPP11 Equilibrium Storage Problem
disp('DEMAPP11 Equilibrium Storage Problem')
close all

% Enter model parameters
r   = 0.1;
k   = 0.5;
eta = 5;
s0  = 1;

% Define the approximant
T = 1;
n = 15;
tnodes = chebnode(n-1,0,T);
fspace = fundefn('cheb',n,0,T);

% Call rootfinding algorithm
c = zeros(n,2); c(1,:) = 1;
c = broyden('fapp11',c(:),tnodes,T,n,fspace,r,k,eta,s0);

% Compute solution and residual functions
nplot = 501;
t = nodeunif(nplot,0,T);
c = reshape(c,n,2);
x = funeval(c,fspace,t);
d = funeval(c,fspace,t,1);
r = d - [r*x(:,1)+k  -x(:,1).^(-eta)];

figure(1)
plot(t,x(:,1),t,x(:,2))
xlabel('Time')
% ylabel('x')
ylim([0 1.6])
legend('x_1(t): Price','x_2(t): Stocks',2)
title('Solution Functions for Equilibrium Storage Problem')

figure(2)
plot(t,r(:,1),t,r(:,2))
xlabel('Time')
ylabel('Residual')
legend('R_1(t)','R_2(t)',4)
title('Residual Functions for Equilibrium Storage Problem')

disp('Maximum Absolute Residuals')
disp(max(abs(r)))

prtfigs(mfilename,'Residual Functions for Equilibrium Storage Problem',2)
prtfigs(mfilename,'Solution Functions for Equilibrium Storage Problem',1)