% BVPRES Residual function used by BVPSOLVE
% See also BVPSOLVE

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function r=bvpres(c,func,params,fspace,tnode,Phi,Phi1,tb,phi,phi1)
  n=size(Phi,2);
  m=length(tb);
  c=reshape(c,n,m);

% Compute residuals at nodal values
  x=Phi*c;
  dx=Phi1*c;
  r=feval(func,'r',tnode,x,dx,params{:});

% Compute boundary conditions and concatenate to residuals
  x=phi*c;
  dx=phi1*c;
  b=feval(func,'b',tb,x,dx,params{:});
  r=[r(:);b(:)];
