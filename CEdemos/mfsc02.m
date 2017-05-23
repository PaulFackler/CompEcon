% MFSC03 Model function file for social planner's renewable resource problem
function [out1,out2]=mfsc03(flag,s,x,Vs,alpha,beta,K,b,eta,C,gamma,sigma,rho)

switch flag
case 'foc'
  Cost=C*s.^(-gamma);
  out1=b*x.^(-1/eta)-Cost-Vs;
  out2=-(b/eta)*x.^(-1/eta-1);
case 'x'
  Cost=C*s.^(-gamma);
  out1=b*(Cost+Vs).^(-eta);
case 'f'
    Cost=C*s.^(-gamma);
    if eta~=1                          % demand elasticity <> 1
      factor1=1-1/eta;              % case separately to avoid 
      factor0=b.^(1/eta)/factor1;  % division by 0; see iteration loop
      out1=factor0*x.^factor1-Cost.*x;
    else                               % demand elasticity = 1
      out1=b*log(x)-Cost.*x;
    end  
case 'g'
  if beta~=0                               % need to handle beta=0
    Growth=alpha/beta*s.*(1-(s/K).^beta);  % case separately
  else                                     % to avoid division by 0
    Growth=alpha*s.*log(K./s);              
  end
  out1=Growth-x;
case 'sigma'
  out1=sigma*s;
case 'rho'
  out1=rho+zeros(size(s,1),1);
end