% BVPRES Used by BVPSOLVE
function r=bvpres(c,modelfunc,modelparams,fspace,tnode,Phi,Phi1,tb,phi,phi1)

  n=size(Phi,2);
  m=length(tb);
  c=reshape(c,n,m);

  % Compute residuals at nodal values
  x=Phi*c;
  dx=Phi1*c;
  r=feval(modelfunc,'r',tnode,x,dx,modelparams{:});

  % Compute boundary conditions and concatenate to residuals
  x=phi*c;
  dx=phi1*c;
  b=feval(modelfunc,'b',tb,x,dx,modelparams{:});
  r=[r(:);b(:)];
