% DEMRS03 Entry/exit Model
% ref.: Dixit & Pindyck, p. 218

function demdsc03 
  close all
  optset('rssolve','defaults')
  
  disp(' ')
  disp('DEMRS03 Entry/exit Model')

% Define parameters 
  r     = 0.05;
  mu    = 0; 
  sigma = 0.2; 
  C     = 1; 
  E     = 2;
  I     = 5; 

% Pack model structure
  clear model
  model.func='mfrs03';
  model.params={r,mu,sigma,C,E,I};
  model.xindex=...
   [1 0 0 0 0;
    1 2 1 2 0;
    2 1 1 2 0;
    2 0 0 0 1];

% Set initial values and approximation size  
  x=[0; 3; 1; 50];
  n=[15 85];

% Call solver
  [cv,fspace,x]=rssolve(model,x,n);

  x=reshape(x,2,2);
  pstar=[x(2,1);x(1,2)];

% Compute approximation value functions 
  P=linspace(0,3,501)';

  V=[funeval(cv{1},fspace{1},P),...
     funeval(cv{2},fspace{2},P)];
  V(P>x(2,1),1)=V(P>x(2,1),2)-I;
  V(P<x(1,2),2)=V(P<x(1,2),1)-E;

  dV=[funeval(cv{1},fspace{1},P,1),...
      funeval(cv{2},fspace{2},P,1)];
  dV(P>x(2,1),1)=dV(P>x(2,1),2);
  dV(P<x(1,2),2)=dV(P<x(1,2),1);
  d2V=[funeval(cv{1},fspace{1},P,2),...
       funeval(cv{2},fspace{2},P,2)];
  d2V(P>x(2,1),1)=d2V(P>x(2,1),2);
  d2V(P<x(1,2),2)=d2V(P<x(1,2),1);

  e=r*V-(mu*P(:,[1 1])).*dV-0.5*sigma^2*P(:,[1 1]).^2.*d2V;
  e(:,2)=e(:,2)-(P-C);
  e(P>x(2,1),1)=nan;
  e(P<x(1,2),2)=nan;

% Plot results
  vstar=[funeval(cv{1},fspace{1},pstar(1));
         funeval(cv{2},fspace{2},pstar(2))];
  figure(1)
  plot(P,V);
  hold on; plot(pstar(1),vstar(1),'*',pstar(2),vstar(2),'*'); hold off
  title('Value Functions')
  xlabel('P')
  xlim([0 3])
  ylabel('V(P)')
  legend('Inactive','Active',4)

  dvstar=[funeval(cv{1},fspace{1},pstar(1),1);
         funeval(cv{2},fspace{2},pstar(2),1)];
  figure(2)
  dV(P>pstar(1),1)=NaN;
  dV(P<pstar(2),2)=NaN;
  plot(P,dV);
  hold on; plot(pstar(1),dvstar(1),'*',pstar(2),dvstar(2),'*'); hold off
  title('Marginal Value Functions')
  xlabel('P')
  ylabel('V''(P)')
  xlim([0 3])
  legend('Inactive','Active',4)

  figure(3)
  plot(P,e)
  title('Approximation Residual')
  xlabel('P')
  ylabel('Residual')
  ylim([-8e-6 8e-6])
  legend('Inactive','Active',4)

% Compute "closed form" solutions
  optset('broyden','showiters',0)
  optset('broyden','maxit',100)
  [residual,beta1,beta2,A1,A2,Pl,Ph] = entexgbm([],r,mu,sigma,I,E,C);
  VVi=(A1*P.^beta1);
  VVa=(A2*P.^beta2+P./(r-mu)-C./r);
  VVi(P>Ph)=VVa(P>Ph)-I;
  VVa(P<Pl)=VVi(P<Pl)-E;

  disp('Boundary points: "Exact", approximate and error')
  disp([[Ph;Pl] pstar [Ph;Pl]-pstar])

prtfigs(mfilename,'Solution to the Entry/Exit Model',[1 2 3])
