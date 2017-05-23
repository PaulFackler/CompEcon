% FAPP11 Residual function for equilibrium storage problem
% See DEMAPP11
function r = fapp11(c,tnodes,T,n,fspace,r,k,eta,s0);
c  = reshape(c,n,2);
x  = funeval(c,fspace,tnodes);
d  = funeval(c,fspace,tnodes,1);
r  = d - [r*x(:,1)+k  -x(:,1).^(-eta)];
x0 = funeval(c,fspace,0);
x1 = funeval(c,fspace,T);
r  = [r(:); x0(2)-s0 ; x1(2)-0];