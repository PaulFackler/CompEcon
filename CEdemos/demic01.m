% DEMIC01 Asset Replacement Model
  disp(' ')
  disp(' DEMIC01 Asset Replacement Model')
  close all

% Define parameters
  beta = [1; 0.05; -0.003];
  P    = 2;
  C    = 3;
  rho  = 0.1;

% Create model variable
  model.func='mfic01';
  model.params={beta,P,rho};
  model.xindex=[0 0;2 0];
  model.F=[0;C];

% Define starting values 
  x=[0 0;100 0];

% Call solver
  n=15;
  optset('broyden','showiters',1)
  [cv,fspace,x]=icsolve(model,x,n);
  Astar=x(2,1);

% Plot results
  A=linspace(0,Astar,101)';
  V=funeval(cv,fspace,A);
  figure(1)
  plot(A,V,'k',Astar,V(end),'k*');
  title('Value Function')
  xlabel('A')
  ylabel('V');
  h=text(15,19,['A^* =' num2str(Astar)]);
  set(h,'HorizontalAlignment','right')
  xlim([0 18])

  prtfigs(mfilename,'Value Function for the Asset Replacement Model',[1])