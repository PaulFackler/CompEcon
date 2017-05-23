% FSC03 Function definition file for neoclassical optimal growth problem
% See DEMSC03 for usage
  function out=fsc03(flag,in,alpha,delta,rho,gam)
  switch flag
  case 'q'
    out=alpha*log(in+1)-delta*in;
  case 's'
    out=alpha./(in+1)-(delta+rho);
  case 'r'
    out=gam./in;
  end