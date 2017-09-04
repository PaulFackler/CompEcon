function [y,p] = qnwlognalt(ny,mu,var)
c = nodeunif(ny+1,0,1);
c = (c(1:ny)+c(2:ny+1))/2;
y = icdf('Lognormal',c,mu,sqrt(var));
y = ny*y*exp(mu+var/2)/sum(y);
p = ones(ny,1)/ny;