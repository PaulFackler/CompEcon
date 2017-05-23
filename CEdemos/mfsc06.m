% MFSC06 Model function file for nonrenewable resource model
function out = mfsc06(flag,s,x,Vs,A,alpha,rho);
 n = size(s,1);
 switch flag
 case 'x';                                           % OPTIMAL CONTROL
   out = (s.^(-alpha/(1-alpha)).*Vs./A).^(-1/alpha);
 case 'f';                                           % REWARD FUNCTION
   out = A*x.^(1-alpha);
 case 'g';                                           % DRIFT FUNCTION
   out = -x.*(1-alpha).*(s.^(-alpha/(1-alpha)));
 case 'sigma'                                        % DIFFUSION FUNCTION
   out = [];
 case 'rho'                                          % DISCOUNT FUNCTION
   out = rho+zeros(n,1);
 end