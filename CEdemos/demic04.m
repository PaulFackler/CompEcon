% DEMIC04 Capacity Choice Model
  close all
  disp(' ')
  disp('DEMIC04 Capacity Choice Model')

% Define parameters
  P     = 2;
  C     = 1;
  delta = 0.5;
  rho   = 0.1;

% Create model variable
  model.func='mfic04';
  model.params={P,C,delta,rho};
  model.xindex=[2 0;0 0];
  model.F=[0;0];
 
% Define starting values 
  xmax=20;
  x=[0 0;xmax 0];

% Call solver
  n=25;
  optset('broyden','showiters',1)
  [cv,fspace,x]=icsolve(model,x,n);
  Kstar=x(1,1);

% Plot results
  K=linspace(0,xmax,101)';
  V=funeval(cv,fspace,K);
  Vstar=funeval(cv,fspace,Kstar);
  V(K<Kstar)=Vstar-(Kstar-K(K<Kstar))*C;
  figure(1)
  plot(K,V,'k',Kstar,Vstar,'k*');
  title('Value Function')
  xlabel('K')
  ylabel('V');
  text(3,Vstar,['K^* =' num2str(Kstar)]) 

  disp('       K*     V(K*)')
  disp([Kstar Vstar])

  prtfigs(mfilename,'Value Function for the Capacity Choice Model')