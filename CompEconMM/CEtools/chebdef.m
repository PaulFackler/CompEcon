%% CHEBDEF 
%
%  Defines parameters for Chebychev polynomial functions
%
%  Usage
%    [n,a,b,parms] = chebdef(n,a,b)
%  Input
%    n         : number of basis functions (1 plus the polynomial order)
%    a         : lower bound of approximation interval
%    b         : upper bound of approximation interval
%  Output
%    n         : number of basis functions (1 plus the polynomial order)
%    a         : lower bound of approximation interval
%    b         : upper bound of approximation interval
%    parms     : cell array containing n, a, b
%  See Also    : chebase

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu

function [n,a,b,parms] = chebdef(n,a,b)

if nargin<3, error('3 parameters must be specified'), end
if n<=0 | fix(n)~=n,
  error('n must be a positive integer')
end
if (a>=b)
  error('Lower bound a must be less than upper bound b');
end
parms = {n,a,b};