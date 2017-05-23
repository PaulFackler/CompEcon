% BVPSOLVE Solves general first order boundary value problems
%    f(t,x(t),x'(t)) = 0
% s.t.
%    b(tb,x(tb),x'(tb)) = 0.
% USAGE
%   [c,x,r] = bvpsolve(model,fspace,tnode,c,t)
% INPUT
%   model  : model structure variable (described below)
%   fspace : function space (defined with FUNDEF)
%   tnode  : n-1 nodal points for collocation
%   c      : initial coefficient values (n by d)
%   t      : (m by 1) vector of evaluation points (optional)
% OUTPUT
%   c      : coefficients satisfying the collocation conditions
%   x      : the solution evaluated at t (m by d)
%   r      : the residual function evaluated at t (m by d)
%
% The model structure has three fields:
%   func   : the name of a model definition file (described below)
%   tb     : a set of d time values at which boundary conditions are
%             evaluated
%   params : parameters that are passed to the func file
%
% The function definition file must have the syntax
%     function out = bvpfile(flag,t,x,dx,additional parameters)
%     switch flag
%     case 'r'
%        out = f(t,x,dx);   % the residual function 
%     case 'b'
%        out = b(tb,x,dx);  % the boundary condition 
%     end
% The additional parameters used by this function should match
%   the set in model.params
%
% USES: BROYDEN

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [c,x,r]=bvpsolve(model,fspace,tnode,c,t)

d=size(c,2);
tb=model.tb;

if length(tb)~=d
  error('The # of boundary conditions must equal the dimension of the system')
end
if fspace.d~=1
  error('BVPSOLVE only solves ODEs')
end
if size(tnode)~=[fspace.n-1,1]
  error('tnode should be n-1 by 1')
end

% nodal basis matrices
  Phi=funbas(fspace,tnode);
  Phi1=funbas(fspace,tnode,1);

% boundary point basis matrices
  phi=funbas(fspace,tb);
  phi1=funbas(fspace,tb,1);

% Call rootfinding algorithm
  c=broyden('bvpres',c(:),model.func,model.params,fspace,tnode,Phi,Phi1,tb,phi,phi1);
  c=reshape(c,fspace.n,d);

% compute solution and residual functions
  if nargout>1 & nargin>4 & ~isempty(t)
    x=funeval(c,fspace,t);
    dx=funeval(c,fspace,t,1);
    r=feval(model.func,'r',t,x,dx,model.params{:});
  end
