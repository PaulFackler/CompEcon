% DEMAPP09 Cournot oligopolist problem
disp('DEMPAPP09 Cournot oligopolist problem')
close all

alpha = 1.0; 
eta   = 1.5;

n =  25; 
a = 0.1;  
b = 3.0;
fspace = fundefn('cheb',n,a,b);
p      = funnode(fspace);

c = ones(n,1);
optset('broyden','defaults')
c = broyden('cournotres',c,p,fspace,alpha,eta);
c = real(c);

figure(1)
pplot = nodeunif(501,a,b);
splot = funeval(c,fspace,pplot);
plot(splot,pplot);
title('Cournot Effective Firm Supply Function')
xlabel('Quantity'); ylabel('Price'); 

figure(2);
rplot = cournotres(c,pplot,fspace,alpha,eta);
plot(pplot,rplot)
title('Residual Function for Cournot Problem')
xlabel('Price'); ylabel('Residual')

figure(3)
m = [1 3 5 10 15 20];

hp=plot(splot*m,pplot,pplot.^(-eta),pplot);
title('Industry Supply and Demand Functions')
xlabel('Quantity'); ylabel('Price')
hl = legend('m=1','m=3','m=5','m=10','m=15','m=20');
xlim([0 13]);

figure(4);
pp = (b+a)/2;
dp = (b-a)/2;
m  = (1:25)';
for i=1:50
  dp = dp/2;
  pp = pp-sign(funeval(c,fspace,pp).*m-pp.^(-eta)).*dp;
end
plot(m,pp)
title('Cournot Equilibrium Price as Function of Industry Size')
xlabel('Number of Firms'); ylabel('Price')

prtfigs(mfilename,'Approximation Residual for the Cournot Model',2)
prtfigs(mfilename,'Cournot Industry Supply and Demand Functions',3)
prtfigs(mfilename,'Cournot Equilibrium Price as Function of Industry Size',4)
