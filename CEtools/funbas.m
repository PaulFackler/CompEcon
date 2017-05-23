% FUNBAS Computes a basis matrix
% USAGE
%   B=funbas(fspace,x,order);
% INPUTS
%   fspace  : a structure defining a family of functions (use FUNDEF to create this)
%   x       : an mxd matrix or a 1xd cell array of columns vectors
%                   (default: created by calling FUNNODE)
%   order   : a 1xd vector (default: zeros(1,d)))
% OUTPUT
%   B       : an m by prod(fspace.n) basis matrix

%
% Note: "order" should consist of a single row.
%
% USES: FUNBASX
%
% See also: FUNBASX

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function B=funbas(fspace,x,order)

if nargin<1 | ~isstruct(fspace)
  error('Must pass a function definition structure');
end
if nargin<3 | isempty(order)
  order=zeros(1,fspace.d);
end
if nargin<2 | isempty(x)
  x=funnode(fspace);
end

if size(order,1)>1,
  warning('In FUNBAS "order" should have only one row, other rows are ignored')
  order=order(1,:);
end

B=funbasx(fspace,x,order,'expanded');
B=B.vals{1};
