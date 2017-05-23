% DIGAMMA Computes the digamma (psi) function for positive arguments
% The digamma (psi) function is the first derivative of the loggamma function.
% USAGE
%   y=digamma(x);
% Accepts positive matrices.
% Based on formula 6.3.18 with recurrence formula 6.3.5 in Abromowitz and Stegun.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function y=digamma(x)
nind=x<=0;
uind=x<1 & ~nind;
xu=1-x(uind);
x(uind)=xu;
x=x+10;
y=1./(x.*x);
y=((((((-8.33333333333333333333e-2*y   ...
        +2.10927960927960927961e-2).*y ...
        -7.57575757575757575758e-3).*y ...
        +4.16666666666666666667e-3).*y ...
        -3.96825396825396825397e-3).*y ...
        +8.33333333333333333333e-3).*y ...
        -8.33333333333333333333e-2).*y;
y=y+log(x)-0.5./x-1./(x-10)-1./(x-9)-1./(x-8)-1./(x-7)-1./(x-6) ...
                 -1./(x-5)-1./(x-4)-1./(x-3)-1./(x-2)-1./(x-1);
y(nind)=NaN;
y(uind)=y(uind)+pi*cot(pi*xu);
