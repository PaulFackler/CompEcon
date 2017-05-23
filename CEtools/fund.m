% FUND Evaluates functions and first 2 derivatives
% USAGE
%   [F,J,H]=fund(c,fspace,x,HessOpt);
% INPUTS
%   c       : a coefficient matrix (n x p)
%   fspace  : a family definition structure
%   x       : a set of evaluate points (m x d)
%   HessOpt : 0/1/2 (default=0; see below)
% OUTPUTS
%   F : the function values
%   J : the Jacobian
%   H : the Hessian
%         
% The Jacobian, J, is an (m x p x d) array where
% J(i,j,k) is the derivative of the jth function with
% respect to the kth input variable evaluated at x(i,:)
%
% The Hessian, H, can take 3 forms depending on HessOpt:
%  HessOpt = 0: vech form  - an (m x p x d(d+1)/2) array 
%          = 1: vec form   - an (m x p x d^2) array 
%          = 2: d x d form - an (m x p x d x d) array
% In general H(i,j,:) is the second derivative of the jth function
% with respect to the input variables evaluated at x(i,:)
%
% See also: FUNJAC, FUNHESS

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [F,J,H]=fund(c,fspace,x,HessOpt)

d=fspace.d;
order=(0:nargout-1)';       % Only evaluate the bases that are needed
B=funbasx(fspace,x,order);

% Evaluate the function
F=funeval(c,fspace,B);
if nargout>1
  % Evaluate the Jacobian
  order=eye(d);
  J=funeval(c,fspace,B,order);
  if nargout>2
    % Evaluate the Hessian
    % Create the order matrix for the Hessian
    u=ones(d,1); eyed=eye(d);
    order=kron(eyed,u)+kron(u,eyed);
    ind=find(tril(ones(d,d))); 
    order=order(ind,:);
    H=funeval(c,fspace,B,order);
    % Expand the Hessian is vech form is not wanted
    if HessOpt>0
      H=H(:,:,vechinv(1:d*(d+1)/2));
      if HessOpt>1
        H=reshape(H,size(H,1),size(H,2),d,d);
      end
    end
  end
end