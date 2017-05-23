% ARRAYSS Reformulates an MCP as a semismooth function
% USAGE
%   [fxnew,Jnew]=arrayss(fx,x,a,b,J);
% INPUTS
%    x : an evaluation point
%    a : lower bounds
%    b : upper bounds
%   fx : function values at x
%    J : Jacobian of f at x (optional, required if Jnew requested)
% OUTPUTS
%   fxnew : value of semi-smooth function at x
%    Jnew : Jacobian of function at x
%
% The reformulation uses
%   phi^-(phi^+(fx,a-x),b-x)
% where
%   phi^+(y,z)=y+z+sqrt(y.^2+z^2)
%   phi^-(y,z)=y+z-sqrt(y.^2+z^2)

% Copyright 1997-2001 by Paul L. Fackler (paul_fackler@ncsu.edu)

function [fxnew,Jnew]=arrayss(x,a,b,fx,J)
if nargout<2
  fxnew=arrayssx(x,a,b,fx);  
else                            % Compute the Jacobian
  [fxnew,ff,aa]=arrayssx(x,a,b,fx);  
  [m,n]=size(x);
  Jnew=repmat(ff,1,n).*reshape(J,m,n*n);
  % index for diagonal elements of the Jacobian
  ind=(0:n:n*(n-1)) + (1:n);
  Jnew(:,ind)=Jnew(:,ind)-aa;
end