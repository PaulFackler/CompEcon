function [out1,out2] = mfre02(flag,s,x,z,e,snext,xnext,delta,alpha,rhoz,rhor,gamma,beta,sigma,omega,bbar);
% Function file for optimal growth model

k=s(:,1);
b=s(:,2);
zvar=exp(s(:,3));
rstar=s(:,4);
c=x(:,1);
n=x(:,2);
i=x(:,3);
r=x(:,4);


switch flag
case 'f'; % REWARD FUNCTION
  uc=gamma*(c.^(gamma-1).*n.^(1-gamma)).^(-sigma);
  un=(1-gamma)*(c.^(gamma).*n.^(-gamma)).^(-sigma);
  mvl=(1-alpha)*exp(zvar).*k.^(alpha).*n.^(-alpha);
  out1 = [uc-beta*z(:,1) ...
          uc-beta*z(:,2) ...
          un-uc.*mvl      ...
          rstar+omega*(b-bbar)-r];
case 'x'
  out1 = [];
case 'g'; % STATE TRANSITION FUNCTION
  out1 = [i+(1-delta)*k ...
          (1+r).*b-i-c+exp(zvar).*k.^alpha.*n.^(1-alpha)...
          rhoz*zvar+e(:,1)...
          rhor*rstar+e(:,2)];           
case 'h'
  kn=snext(:,1);
  bn=snext(:,2);
  zn=snext(:,3);
  rstarn=snext(:,4);
  cn=xnext(:,1);
  nn=xnext(:,2);
  in=xnext(:,3);
  rn=xnext(:,4);

  uchat=gamma*(cn.^(gamma-1).*nn.^(1-gamma)).^(-sigma);
  mvk=alpha*exp(zn).*kn.^(alpha-1).*nn.^(1-alpha);
  out1= [uchat.*(1+rstarn) uchat.*mvk];
end 