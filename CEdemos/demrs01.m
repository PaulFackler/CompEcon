% DEMRS01 Asset Abandonment Model
  close all
  optset('rssolve','defaults')

  disp(' ')
  disp('DEMRS01 Asset Abandonment Model')
  
% Define parameters
  c     = 0.5;
  mu    = 0;  
  sigma = 0.2;  
  rho   = 0.1;  
  
% Pack model structure
  clear model
  model.func='mfrs01';
  model.params={c,mu,sigma,rho};
  model.xindex=[1 0 1 2 0;1 0 0 0 1];

% Set approximation size and initial values
  n=101;
  x=[0;20];
  
% Call solver
  [cv,fspace,x]=rssolve(model,x,n,'cheb');
  cv=cv{1};
  fspace=fspace{1};
  
% Plot results
  S=linspace(0,1,1001)';
  V=funeval(cv,fspace,S);
  Vs=funeval(cv,fspace,S,1);
  Vss=funeval(cv,fspace,S,2);

  sstar=x(1);

  V(S<sstar)=0;  
  Vs(S<sstar)=0;  
  Vss(S<sstar)=0;

  figure(1)
  plot(S,V,'k',sstar,funeval(cv,fspace,sstar),'k*');
  title('Value Function')
  xlabel('P')
  ylabel('V');

  figure(2)
  Vs=funeval(cv,fspace,S,1);
  Vs(S<sstar)=0;
  plot(S,Vs,'k',sstar,funeval(cv,fspace,sstar,1),'k*');
  title('Marginal Value Function')
  xlabel('P')
  ylabel('V''');

  S=linspace(0,fspace.b,1001)';
  V=funeval(cv,fspace,S);
  Vs=funeval(cv,fspace,S,1);
  Vss=funeval(cv,fspace,S,2);

  sstar=x(1);

  V(S<sstar)=0;  
  Vs(S<sstar)=0;  
  Vss(S<sstar)=0;

% Exact solution
  beta=roots([sigma^2/2 mu-sigma^2/2 -rho]);
  beta=beta(beta<0);
  pstar=(rho-mu)*beta*c/rho/(beta-1);
  A=-pstar.^(1-beta)/(rho-mu)/beta;
  VV=S/(rho-mu)-c/rho+A*S.^beta;
  VV(S<pstar)=0;

  e=rho*V-(S-c).*(S>sstar)-mu*S.*Vs-0.5*sigma^2*S.^2.*Vss;
 
  figure(3) 
  plot(S,e);
  title('Aproximation Residual')
  xlabel('P')
  ylabel('Residual')

  figure(4)
  plot(S,VV-V(:,end));
  title('Approximation Error')
  xlabel('P')
  ylabel('Error')

  disp('P*: Approx    Exact     Error')
  disp([sstar pstar sstar-pstar])

prtfigs('demrs01','Solution to the Asset Abandonment Model',[1 2 3 4])