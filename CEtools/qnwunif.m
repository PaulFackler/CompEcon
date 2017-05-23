% QNWUNIF Computes nodes and weights for multivariate uniform distribution
% USAGE
%   [x,w] = qnwunif(n,a,b);
% INPUTS
%   n   : 1 by d vector of number of nodes for each variable
%   a   : 1 by d vector of left endpoints for interval
%   b   : 1 by d vector of right endpoints for interval
% OUTPUTS
%   x   : prod(n) by d matrix of evaluation nodes
%   w   : prod(n) by 1 vector of probabilities
% 
% To compute expectation of f(x), where x is U(a,b), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwunif(n,a,b); E[f]=w'*f(x);
%
% USES: qnwlege.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu 

function [x,w] = qnwunif(n,a,b);

[x,w] = qnwlege(n,a,b);
w = w/prod(b-a);
