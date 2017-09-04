%% LPMETRIC
%
%  Computes Lp distance between two real-valued functions defined on bounded interval
%
%
%  Usage
%    y = lpmetric(f,g,p,l,u)
%  Input
%    f         : real-valued function of form y=f(x)
%    g         : real-valued function of form y=f(x)
%    p         : p>0 or inf, order of metric
%    l         : lower bound of integration interval
%    u         : upper bound of integration interval
%  Output
%    y         : Lp distance between f and g

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function y = lpmetric(f,g,p,l,u)

if p<=0, error('p must be positive or infinite'), end

if p<inf
    y = integral(@(x)abs(f(x)-g(x)).^p,l,u)^(1/p);
else
    x = nodeunif(500,l,u);
    y = max(abs(f(x)-g(x)));
end