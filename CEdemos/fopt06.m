function [fx,fjac] = fopt06(x)
fx   = x.*cos(x.^2);
fjac = cos(x.^2)-2*x.^2.*sin(x.^2);
