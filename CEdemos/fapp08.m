% FAPP08 A 3-D function to be approximated
function [y,d]=fapp08(x)
y = (exp(x(:,1).*x(:,2))+x(:,3).^3)/100;
d(:,1) = (x(:,2).*exp(x(:,1).*x(:,2)))/100;
d(:,2) = (x(:,1).*exp(x(:,1).*x(:,2)))/100;
d(:,3) = (3*x(:,3).^2)/100;
