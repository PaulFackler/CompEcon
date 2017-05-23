function [fval,fjac] = billups(x)
% Billups' function
fval = 1.01-(1-x).^2;
fjac = 2*(1-x);
