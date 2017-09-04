%% FUNNODE
%
%  Generates default interpolation nodes for a family of basis functions
%
%  Usage
%    [x,xcoord] = funnode(basis)
%  Let
%    d = dimension of basis interval domain
%    n = 1.d number of coordinate nodes per dimension
%    a = 1.d interval lower bounds
%    b = 1.d interval upper bounds
%  Input
%    basis     : d-dimensional family of basis functions created with fundefn
%  Output
%    x         : prod(n).d grid of interpolation nodes on [a,b]
%    xcoord    : 1.d cell array of node coordinates by dimension
%  Note
%    If d=1, x is an by by 1 vector of nodes on the interval [a,b] in R.
%    If d>1, x is an N=prod(n) by d grid of nodes on the interval [a,b] in
%    R^d and xcoord is a by 1 cell array where xccord{i} is a n(i) by 1
%    vector of coordinates along dimension i.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,xcoord] = funnode(basis)

d = basis.d;
if d==1
  xcoord = feval([basis.bastype{1} 'node'],basis.parms{1}{:});
  x = xcoord;
else
  xcoord = cell(1,d);
  for j=1:d
    xcoord{j} = feval([basis.bastype{j} 'node'],basis.parms{j}{:});
  end
  x = gridmake(xcoord);
end