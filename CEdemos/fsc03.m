% FSC03 Residual function for neoclassical optimal growth problem
% See DEMSC03 for usage
function e=fsc03(c,r,svals,qvals,Phi0,Phi1,Cstar)
C=Phi0*c;
dC=Phi1*c;
warning off
e=(qvals-C).*dC-svals./feval(r,C);
warning on
e(1)=C(1);
e(2)=C(2)-Cstar;