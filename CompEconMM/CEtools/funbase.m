%% FUNBASE
%
%  Computes basis function values and derivatives at prescribed points.
%
%  Usage
%    B = funbase(basis,x,order)
%  Let
%    m = number of evaluation points
%    d = dimension of evaluation points
%    N = prod(basis.n) number of basis functions
%  Input
%    basis  : family of d-dimensional basis functions created with fundefn
%    x      : m.d matrix or 1xd cell array of m.1 column vectors
%    order  : 1.d vector indicating order of differentiation
%  Output
%    B      : m.N basis matrix
%  Note
%    If x not passed, x set to N.d matrix of standard basis nodes.
%    If order not passed, order set to zeros(1,d) and function values are computed.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function B = funbase(basis,x,order)

if nargin<1 || ~isstruct(basis)
  error('In funbase: "basis" must be a structure');
end
if nargin<3 || isempty(order)
  order = zeros(1,basis.d);
end
if nargin<2 || isempty(x)
  x = funnode(basis);
end
if size(order,1)>1,
  warning('In funbase: "order" should have only one row')
  order = order(1,:);
end
B = funbasex(basis,x,order,'expanded');
B = B.vals{1};