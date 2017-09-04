%% QNWLOGN
%
%  Discretizes lognormal distribution with log mean mu and log variance var.
%
%  Equivalently, discretizes random variable whose natural logarithm is
%  normally distributed with mean mu and variance var.
%
%  Usage
%    [x,w] = qnwlogn(n,mu,var)
%  Input
%    n   : number of nodes per dimension
%    mu  : log mean (default: zeros)
%    var : positive log variance (default: 1)
%  Output
%    x   : n.1 discrete mass points
%    w   : n.1 associated probabilities
%  Note
%    To compute Ef(X) when f is real-valued and X is Lognormal(mu,var), 
%    write a Matlab function f that returns n.1 vector when passed an m.1 
%    vector, and execute [x,w]=qnwlogn(n,mu,var); Ef=w'*f(x).
%  Note
%    The lognormal distribution is defined on (0,inf)^d.  The mean and
%    and variance of the univariate distribution are exp(mu+var/2) and 
%    (exp(var)-1)exp(2mu+var), respectively.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwlogn(n,mu,var)

if nargin<2, mu  = 0; end
if nargin<3, var = 1; end

[x,w] = qnwnorm(n,mu,var);
x = exp(x);