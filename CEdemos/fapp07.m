% FAPP07 A 2-D function to be approximated
function [y,d]=fapp07(x)
y=x(:,1)./exp(x(:,2));
d(:,1)=1./exp(x(:,2));
d(:,2)=-x(:,1)./exp(x(:,2));
