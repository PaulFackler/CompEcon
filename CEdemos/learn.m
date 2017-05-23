% Learn Solves sequential learning problem
function [Q,pstar,c,cdef,A1,A2,beta1,beta2]=learn(r,delta,sigma,Qm,cbar,C,N,n,b)
% Compute solution for Q>Qm
  beta1=0.5-(r-delta)/sigma.^2 + sqrt((0.5-(r-delta)/sigma.^2).^2+2*r/sigma.^2);
  beta2=0.5-(r-delta)/sigma.^2 - sqrt((0.5-(r-delta)/sigma.^2).^2+2*r/sigma.^2);
% Impose value matching and smooth pasting at P=cbar
  temp=[       cbar.^beta1         -(cbar.^beta2)     ; ...
        beta1*cbar.^(beta1-1)   -beta2*cbar.^(beta2-1)];
  temp=temp\[cbar/delta-cbar/r ; 1/delta];
  A1=temp(1);  A2=temp(2);
% Define the approximating functions and nodal values
  Delta=Qm/N;
  Q=linspace(0,Qm,N+1);  
  cdef=fundefn('cheb',n-1,0,b);
  y=funnode(cdef);
  cdef=fundefn('cheb',n,0,b);
% Set up collocation matrices
  D=funbasx(cdef,y,[0;1;2]);
  Phi0=D.vals{1};
  Phi1=D.vals{2};
  D=r*Phi0-(r-delta-0.5*sigma^2)*Phi1-0.5*sigma.^2*D.vals{3};
  phi=funbasx(cdef,0,[0;1;2]);
  A=[           Phi0               zeros(n-1,1);
     phi.vals{2}-beta1*phi.vals{1}       0       ;
     phi.vals{3}-(beta1)*phi.vals{2}     0       ];
  B=[Phi0-Phi1-Delta*D  Delta*exp(y);
              zeros(2,n+1)          ];

  d=[ones(n-1,1);0;0];
% Compute cost function values
  gamma=log(C/cbar)/Qm;
  Cost=Delta*cbar*exp(gamma*(Qm-Q));
% Initialize at terminal boundary
  p=[cbar;cbar*exp(y)];
  c=zeros(n+1,N+1);
  c(1:n,N+1)=[phi.vals{1};Phi0]\(A2*p.^beta2+p/delta-cbar/r);
  c(end,N+1)=cbar;
% Iterate backwards in Q
  for i=N+1:-1:2
    A(1:n-1,end)=Phi1*c(1:n,i)/(-c(end,i));
    c(:,i-1)=A\(B*c(:,i)-Cost(i)*d);
  end
% Extract Pstar
  pstar=c(end,:)';
  c(end,:)=[];