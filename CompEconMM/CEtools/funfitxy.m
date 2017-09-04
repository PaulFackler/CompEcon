%% FUNFITXY
%
%  Computes basis function coefficients for approximant fit to data x and y. 
%
%  Usage
%    c = funfitxy(basis,x,y)
%  Let
%    dx = dimension of x (domain)
%    dy = dimension of y (range)
%    m  = number of data points
%    N  = number of basis functions, N = prod(n) where n(i) is the number
%         of basis functions along dimension i.
%  Input
%    basis     : basis structure created with fundefn
%    x         : m.dx matrix
%    y         : m.dy matrix
%  Output
%    c         : N.dx coefficient matrix
%  Note
%    Once c has been computed, the approximant can be evaluated at any X
%    using funeval(c,basis,X).
%  Note
%    The number of data points m must be the same for x and y and must be
%    equal to or greater than the number of basis functions N.  If N=m, the
%    approximant exactly interpoolstes the data points; if N>m, the
%    coefficients of the least squares approximant are computed.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function c = funfitxy(basis,x,y)

if nargin~=3, error('Three parameters must be specified'); end
m = size(y,1);
if size(x,1)~=m
  error('funfitxy: x and y are incompatible')
end
if prod(basis.n)>m
  error('funfitxy: insufficient number of data points')
end
c = funbase(basis,x)\y;