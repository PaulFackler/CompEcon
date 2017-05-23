function [fx,J]=simple(x);
fx=x-sqrt(x);
J=1-0.5./sqrt(x);
