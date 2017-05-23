% COURNOT Evaluates the equilibrium condition in the Cournot model
function [fval,fjac] = cournot(q)
 c = [0.6; 0.8];
 eta = 1.6; 
 e = -1/eta;
 fval = sum(q)^e + e*sum(q)^(e-1)*q - diag(c)*q;
 if nargout>1
   fjac = e*sum(q)^(e-1)*ones(2,2) + e*sum(q)^(e-1)*eye(2) ...
           + (e-1)*e*sum(q)^(e-2)*q*[1 1] - diag(c);
 end
