% MFSC03 Model function file for Production-Adjustment Model
function out1 = mfsc03(flag,s,x,Vs,kappa,alpha,sigma,a,c,rho);
 n = size(s,1);
 switch flag
 case 'x';                          % OPTIMAL CONTROL
   out1 = Vs(:,2)/a;
 case 'f';                          % REWARD FUNCTION
   out1 = s(:,1).*s(:,2) - 0.5*c*s(:,2).^2 - 0.5*a*x.^2;
 case 'g';                          % DRIFT FUNCTION
   out1 = [kappa*(alpha-s(:,1)) x];
 case 'sigma'                       % DIFFUSION FUNCTION
   out1 = [sigma*s(:,1) zeros(n,3)];
 case 'rho'                         % DISCOUNT FUNCTION
   out1 = rho+zeros(n,1);
 end