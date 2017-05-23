% FUNCONV Converts from one basis family to another
% USAGE
%   c1=funconv(c0,cdef0,cdef1,order)
% INPUTS
%   c0        : a coefficient matrix
%   cdef0     : a family definiton structure (use FUNDEF to create this)
%   cdef1     : new family definition structure
%   order     : order of differentiation (1xd vector) (default: 0)
%  OUTPUT
%    c1       : new coefficient matrix
%
% The approximating function can then be evaluated using FUNEVAL(c1,cdef1,x)
%
% USES: FUNNODE, FUNBASX, CKRONXI
%
% See also: FUNFITXY, FUNDEF, FUNEVAL, FUNNODE, FUNBAS.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function c1=funconv(c0,cdef0,cdef1,order)
if nargin<3, error('Three parameters must be specified'); end
if nargin<4, order=zeros(1,cdef0.d); end 

if cdef0.d~=cdef1.d
  error('Function families must have the same dimension (d)')
end

x=funnode(cdef1);
B=funbasx(cdef1,x);
c1=ckronxi(B.vals,funeval(c0,cdef0,x,order),cdef1.d:-1:1);
