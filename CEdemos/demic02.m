% DEMIC02 Timber Harvesting Model
  close all
  disp(' ')
  disp('DEMIC02 Timber Harvesting Model')

% Define parameters
  alpha = 0.1;
  m     = 1;
  sigma = 0.05;
  P     = 3;
  C     = 0.15;
  rho   = 0.1;

% Create model variable
  model.func='mfic02';
  model.params={alpha,m,sigma,P,rho};
  model.xindex=[0 0;2 0];
  model.F=[0;C];
 
% Define starting values 
  x=[0 0;0.5 0];

% Call solver
  n=15;
  optset('broyden','showiters',1)
  [cv,fspace,x]=icsolve(model,x,n);
  Sstar=x(2,1);

% Plot results
  S=linspace(0,m/2,101)';
  V=funeval(cv,fspace,S);
  Vstar=funeval(cv,fspace,Sstar);
  dV=funeval(cv,fspace,S,1);
  dVstar=funeval(cv,fspace,Sstar,1);
  V(S>Sstar)=Vstar+P*(S(S>Sstar)-Sstar);
  dV(S>Sstar)=P;
  figure(1)
  plot(S,V,'k',Sstar,Vstar,'k*');
  title('Value Function')
  xlabel('S')
  ylabel('V');
  ylim([1.8 3.2])

  figure(2)
  plot(S,dV,'k',Sstar,dVstar,'k*');
  title('Marginal Value')
  xlabel('S')
  ylabel('V''');
  ylim([1.8 3.2])

  prtfigs(mfilename,'Solution to the Timber Harvesting Model',[1 2])