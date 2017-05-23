function [fx,J]=simple(x);
fx=sqrt(x);
J=0.5./sqrt(x);
