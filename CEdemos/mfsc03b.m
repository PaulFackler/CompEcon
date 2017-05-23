% MFSC03b Supplemental model function file for Production-Adjustment Model
% Used to solve for stationary density of price
function out = mfsc03b(flag,s,kappa,alpha,sigma,a,c,rho);
 switch flag
 case 'mu';                        % STATE TRANSITION FUNCTION
   out = kappa*(alpha-s);
 case 'sigma'                      % DIFFUSION FUNCTION
   out = sigma*s;
 end
