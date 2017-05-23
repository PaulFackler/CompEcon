% FUNHESSB Computes basis matrices to evaluate Hessians
% USAGE:
%   [B,order]=funhessb(gdef,x,UnVech)
% INPUTS:
%   gdef    : a family definition structure
%   x       : a set of evaluate points (m x d)
%   UnVech  : optional 0/1 returns full Hessian if = 1
% OUTPUTS:
%   B       : a basis structure
%   order   : the complete set of order indexes
%             
% The Hessian can be evaluated for coefficient matrix g using 
%   H=funeval(g,B,order);
%
% The resulting Hessian, H, is an (m x p x d(d+1)/2) array:
%   H(i,j,:) is the second derivative of the jth function with
%       respect to the input variables evaluated at x(i,:)
%       arranged in vech form
%  Alternatively H is an (m x p x d^2) array if UnVech=1
% 
% use reshape(H) to remove unwanted dimensions
% use vechinv(squeeze(H)) to create a dxd matrix when H is 1x1xd

% Copyright (c) 1999, 2000 by Paul L. Fackler & Mario J. Miranda

function [B,order]=funhessb(gdef,x,format,UnVech)
if nargin<4 | isempty(UnVech), UnVech=0; end
if nargin<3, format=[]; end

d=gdef.d;
% Create the full order
u=ones(d,1); eyed=eye(d);
order=kron(eyed,u)+kron(u,eyed);

if ~UnVech                   % remove redundant evaluations
  ind=find(tril(ones(d,d))); 
  order=order(ind,:);
end
B=funbasx(gdef,x,order,format);
