% FUNJACB Computes basis structure to evaluate Jacobians
% USAGE:
%   [B,order]=funjacb(gdef,x);
% INPUTS:
%   gdef    : a family definition structure
%   x       : a set of evaluate points (m x d)
% OUTPUTS:
%   B       : a basis structure
%   order   : the complete set of order indexes
%         
% The Jacobian can be evaluated for coefficient matrix g using 
%   J=funeval(g,B,order);
% producing an (m x p x d) array where
% J(i,j,k) is the derivative of the jth function with
% respect to the kth input variable evaluated at x(i,:)
%
% See Also: FUNJAC, FUNHESSB

% Copyright (c) 1999, 2000 by Paul L. Fackler & Mario J. Miranda

function [B,order]=funjacb(gdef,x)
  order=eye(gdef.d);
  B=funbasx(gdef,x,order);