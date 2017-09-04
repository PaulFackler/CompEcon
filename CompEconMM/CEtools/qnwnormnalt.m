function [y,p] = qnwnormnalt(ny,mu,var)
c = nodeunif(ny+1,0,1);
y = icdf('Normal',(c(1:ny)+c(2:ny+1))/2,0,1);
y = mu + sqrt(ny*var/sum(y.^2))*y;
p = ones(ny,1)/ny;