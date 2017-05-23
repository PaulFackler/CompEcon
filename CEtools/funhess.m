% FUNHESS Computes the Hessian for FUN functions 
% USAGE
%   H=funhess(c,fspace,x);
% INPUTS
%   c       : a coefficient matrix (n x p)
%   fspace  : a family definition structure
%   x       : a set of evaluate points (m x d)
% OUTPUTS
%   H       : an (m x p x d(d+1)/2) array
%             H(i,j,:) is the second derivative of the jth function
%             with respect to the input variables evaluated at x(i,:)
%             arranged in vech form
% use squeeze(H) to remove unwanted dimensions
% use vechinv(squeeze(H)) to create a dxd matrix when H is 1x1xd
%
% This function sets up an appropriate ORDER matrix and calls FUNEVAL. 
% For example, if d=2, FUNHESS performs the same operation as
%   funeval(c,fspace,x,[2 0;1 1;0 2])

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function H=funhess(c,fspace,x)

d=fspace.d;
% Create the full order
u=ones(d,1); eyed=eye(d);
order=kron(eyed,u)+kron(u,eyed);
% Remove redundent orders (from symmetry)
ind=find(tril(ones(d,d))); 
order=order(ind,:);
% Call funeval
H=funeval(c,fspace,x,order);