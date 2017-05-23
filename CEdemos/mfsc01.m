% MFSC01 Model function file for optimal growth model
function out = mfsc01(flag,s,x,Vs,alpha,delta,gam,rho);
switch flag
case 'x';                         % OPTIMAL CONTROL
  out = Vs.^(-1./gam);
case 'f';                         % REWARD FUNCTION
  out = (x.^(1-gam)-1)./(1-gam);
case 'g';                         % DRIFT FUNCTION
  out = alpha*log(s+1)-delta*s-x;
case 'sigma'                      % DIFFUSION FUNCTION
  out = [];
case 'rho'                        % DISCOUNT FUNCTION
  out = rho+zeros(size(s,1),1);
end