% FAPP09 Residual function for Cournot oligopoly example
function resid=fapp09(c,p,alpha,eta,B);
dp    = (-1./eta)*p.^(eta+1);
q     = max(0,B*c);               % avoid negative quantities
MR    = p+q.*dp;  
MC    = alpha*sqrt(q) + q.^2;
resid = MR-MC;
