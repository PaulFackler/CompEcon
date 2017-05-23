% DEMIC03 Storage Management Model
  close all
  disp(' ')
  disp('DEMIC03 Storage Management Model')

% Define parameters
  mu    = -0.2;
  sigma = 0.05;
  k     = 0.05;
  P     = 2;
  F     = 3;
  rho   = 0.1;

% Create model variable
  model.func='mfic03';
  model.params={mu,sigma,k,P,rho};
  model.xindex=[1 1;0 0];
  model.F=[F;0];
 
% Define starting values 
  smax=8;
  x=[0 smax;smax smax];

% Call solver
  n=10;
  optset('broyden','showiters',1)
  [cv,fspace,x]=icsolve(model,x,n);
  Sstar=x(1,2);

% Plot results
  S=linspace(0,smax,101)';
  V=funeval(cv,fspace,S);
  Vstar=funeval(cv,fspace,Sstar);
  figure(1)
  plot(S,V,'k',Sstar,Vstar,'k*',S,Vstar+P*(S-Sstar));
  title('Value Function')
  xlabel('S')
  ylabel('V');
  text(12,19,['S^* =' num2str(Sstar)]) 
  
  disp('     S         V')
  disp([0 V(1);Sstar Vstar])

  prtfigs(mfilename,'Value Function for the Storage Management Model',[1])