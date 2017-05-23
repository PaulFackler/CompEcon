% QNWLOGN Computes Gauss-Hermite nodes and weights multivariate lognormal distribution
% USAGE
%   [x,w] = qnwlogn(n,mu,var);
% INPUTS
%   n   : 1 by d vector of number of nodes for each variable
%   mu  : 1 by d mean vector
%   var : d by d positive definite covariance matrix
% OUTPUTS
%   x   : prod(n) by d matrix of evaluation nodes
%   w   : prod(n) by 1 vector of probabilities
% 
% To compute expectation of f(x), where log(x) is N(mu,var), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=gqnlogn(n,mu,var); Ef=w'*f(x);
%
% USES: qnwnorm

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,w] = qnwlogn(n,mu,var)

d = length(n);
if nargin<3, var=eye(d); end
if nargin<2, mu=zeros(1,d); end
[x,w] = qnwnorm(n,mu,var);
x = exp(x);
