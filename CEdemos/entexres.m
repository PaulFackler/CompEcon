% ENTEXRES Computes residual function for the entry/exit problem
% See ENTEX

function [r,ca,ci]=entexres(pstar,f,parms,yi,ya,...
    Phii,Phia,phi0i,phi1i,phi1a,phiha,b0,b1,R,cdefi,cdefa);

if isa(parms,'double'), parms={parms}; end

Pl=pstar(1); Ph=pstar(2);
ni=cdefi.n;
na=cdefa.n;

[r,m,s]=feval(f,Ph*yi,parms{:});
Di=DMatrix(Phii,ni,r,m/Ph,s/Ph);
[r,m,s]=feval(f,Pl*ya,parms{:});
Da=DMatrix(Phia,na,r,m/Pl,s/Pl);

B=[[phi0i;Di;phi1i]     zeros(ni,na);
     zeros(na,ni)    [phi1a;Da;phiha]];
B(ni+1,1:ni)=funbas(cdefi,Pl/Ph);
B(ni,ni+1:end)=funbas(cdefa,Ph/Pl);

c=B\(b0+Pl*b1);

R(1,ni+1:end)=Ph/Pl*funbas(cdefa,Ph/Pl,1);
R(2,1:ni)=Pl/Ph*funbas(cdefi,Pl/Ph,1);
r=R*c;

ci=c(1:ni);
ca=c(ni+1:end);

% Compute Differential Operator
function D=DMatrix(Phi,n,r,m,s)
  s=0.5*s.^2;
  u=ones(1,n);
  D=r(:,u).*Phi.vals{1}-m(:,u).*Phi.vals{2}-s(:,u).*Phi.vals{3};
