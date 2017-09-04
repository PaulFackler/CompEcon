%% NODEUNIF
%
%  Generates uniformly spaced nodes for intervals in R^d
%
%  Usage
%    [x,xcoord] = nodeunif(n,a,b)
%  Let
%    d = dimension of coordinate space
%    N = prod(n) number of grid points
%  Input
%    n         : 1.d number of nodes per dimension
%    a         : 1.d interval lower bounds
%    b         : 1.d interval upper bounds
%  Output
%    x         : N.d grid of equally spaced nodes on [a,b]
%    xcoord    : 1.d cell array of node coordinates per dimension
%  Note
%    If d>1, xccord{i} is n(i).1 vector of coordinates along dimension i.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,xcoord] = nodeunif(n,a,b)

d = length(n);
if d==1
  x = linspace(a,b,n)';
  xcoord = x;
else
  for k=1:d
    xcoord{k} = linspace(a(k),b(k),n(k))';
  end
  x = gridmake(xcoord);
end