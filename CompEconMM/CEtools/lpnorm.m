%% LPNORM
%
%  Computes Lp norm of a real-valued function defined on bounded interval
%
%  Usage
%    y = lpnorm(f,p,l,u)
%  Input
%    f         : real-valued function of form y=f(x)
%    p         : p>0 or inf, order of norm
%    l         : lower bound of integration interval
%    u         : upper bound of integration interval
%  Output
%    y         : Lp norm of f

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu

function y = lpnorm(f,p,l,u)

if p<=0, error('p must be positive or infinite'), end

if p<inf
  y = integral(@(x)abs(f(x)).^p,l,u)^(1/p);
else
  x = nodeunif(500,l,u);
  y = max(abs(f(x)));
end