%% DEMAPP09 Monopolist's Effective Supply

% Preliminary tasks
demosetup(mfilename)


%% Residual Function
resid = @(c,basis,p) ...
   p + funeval(c,basis,p)./(-3.5*p.^(-4.5)) ...
     - sqrt(funeval(c,basis,p)) ...
     - funeval(c,basis,p).^2;

% Approximation structure
n = 21; a = 0.5; b = 2.5;
basis = fundefn('cheb',n,a,b);
pnode = funnode(basis);

% Solve for effective supply function
c = [2;zeros(n-1,1)];
c = broyden(resid,c,basis,pnode);

% Setup plot
nplot = 1000;
pplot = nodeunif(nplot,a,b);
rplot = resid(c,basis,pplot);

% Plot effective supply
figure
plot(funeval(c,basis,pplot),pplot)
title('Monopolist''s Effective Supply Curve')
xlabel('Quantity')
ylabel('Price')

% Plot residual
figure
hold on
plot(pplot,rplot) 
plothdash([],0)
title('Collocation Residual')
xlabel('Price')
ylabel('Residual')

%% SAVE FIGURES
printfigures(mfilename)