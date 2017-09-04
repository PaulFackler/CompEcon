%% DEMAPP06 Cournot Oligopoly Model
%
% Solves Cournot oligopoly model for an arbitrary number of m firms using
% Chebychev collocation.

% Preliminary tasks
demosetup(mfilename)

% Model parameters
alpha = 1.0;
eta   = 3.5;
  
% Approximation structure
n =  21;
a = 0.5;
b = 2.0;
basis = fundefn('cheb',n,a,b);
pnode = funnode(basis);

% Residual function
resid = @(c,p) ...
  p + funeval(c,basis,p).*((-1./eta)*p.^(eta+1)) ...
    - alpha*sqrt(funeval(c,basis,p)) ...
    - funeval(c,basis,p).^2;

% Solve for effective supply function
c = [1;zeros(n-1,1)];
c = broyden(resid,c,pnode);

% Plot demand and effective supply for m=5 firms
figure
pplot = nodeunif(501,a,b);
splot = funeval(c,basis,pplot);
dplot = pplot.^-eta;
plot(5*splot,pplot,dplot,pplot);
xlim([0 4])
ylim([0.5 2])
legend('Supply','Demand')
title('Cournot Effective Firm Supply Function')
xlabel('Quantity')
ylabel('Price')

% Plot residual
figure
hold on
rplot = resid(c,pplot);
plot(pplot,rplot)
plothdash([],0)
title('Residual Function for Cournot Problem')
xlabel('Price')
ylabel('Residual')

% Plot demand and effective supply for m=5,10,20 firms
figure
hold on
m = [5 10 20];
plot(splot*m,pplot,pplot.^(-eta),pplot)
legend('Supply m=5','Supply m=10','Supply m=20','Demand')
title('Industry Supply and Demand Functions')
xlabel('Quantity')
ylabel('Price')
xlim([0 13]);

% Compute equilibrium market price as a function of number of firms m
M = 25;
pp = zeros(M,1);
optset('broyden','showiters',0)
for m=1:M
  pp(m) = broyden(@(p) funeval(c,basis,p).*m-p.^(-eta),1);
end

% Plot equilibrium market price as a function of number of firms m
figure
plot((1:M)',pp)
title('Equilibrium Price as Function of Industry Size')
xlabel('Number of Firms');
ylabel('Market Price')


%% SAVE FIGURES
printfigures(mfilename)