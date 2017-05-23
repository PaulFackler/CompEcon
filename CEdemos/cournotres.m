% COURNOTRES Residual function for Cournot oligopoly example
function r = cournotres(c,p,fspace,alpha,eta);
q  = funeval(c,fspace,p);
r  = p+q.*[(-1./eta)*p.^(eta+1)] - alpha*sqrt(q) - q.^2;
