function [y,p] = qnwlognalt2(ny,mu,var)
[y,p] = qnwnormnalt(ny,mu,var);
y = exp(y);
y = ny*y*exp(mu+var/2)/sum(y);