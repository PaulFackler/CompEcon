% DEMIC05 Cash Management Model
  close all
  disp(' ')
  disp('DEMIC05 Cash Management Model')

% Define parameters
  mu    = 0;
  sigma = 0.2;
  rho   = 0.04;
  r     = 0.05;
  c     = 0.75;
  C     = 0.25;
  f     = 1;
  F     = 0.5;

% Create model variable
  model.func='mfic05';
  model.params={mu,sigma,rho,r,c,C};
  model.xindex=[1 1;2 1];
  model.F=[f;F];
 
% Define starting values 
  x=[0 1;2 1];

% Call solver
  n=10;
  optset('broyden','showiters',1)
  [cv,fspace,x]=icsolve(model,x,n);

% Plot results
  S=linspace(0,x(2,1),101)';
  V=funeval(cv,fspace,S);
  Vstar=funeval(cv,fspace,x(:));
  figure(1)
  plot(S,V,'k',x(:),Vstar,'k*');
  title('Value Function')
  xlabel('S')
  ylabel('V');


  dV=funeval(cv,fspace,S,1);
  dVstar=funeval(cv,fspace,x(:),1);
  figure(2)
  plot(S,dV,'k',[0;3],ones(2,1)*(r/rho+[c -C]),'k-',x(:),dVstar,'k*');
  title('Marginal Value Function')
  xlabel('S')
  ylabel('V''')

  % Display selected values
  format short
  disp('         x      V(x)     V''(x)')
  disp([x([1;3;4;2]) Vstar([1;3;4;2]) dVstar([1;3;4;2])])

  %prtfigs(mfilename,'Solution to the Cash Management Model',[1 2])