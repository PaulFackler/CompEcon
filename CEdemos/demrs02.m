% DEMRS02 Optimal Fish Harvest Model
  close all
  optset('rssolve','defaults')

  disp(' ')
  disp('DEMRS02 Optimal Fish Harvest Model')

% Define parameters
  alpha = 0.5; 
  sigma = 0.5;
  H     = 1; 
  P     = 1; 
  c     = 0.25;   
  rho   = 0.1;

% Pack model structure
  clear model
  model.func='mfrs02';
  model.params={alpha,sigma,H,P,c,rho};
  model.xindex=[
    1 0 0 0 0;
    1 2 1 1 2;
    2 0 0 0 0];

% Set initial values and approximation size  
  x=[0.001;0.5;3];
  n=[75 20];

% Call solver
  [cv,fspace,x]=rssolve(model,x,n);

  sstar=x(2);
  aa=[x(1) sstar];
  rr=[sstar x(3)]-aa;

% Plot results
  s=(linspace(0.001,1,2001)'); 
  m=2;
  Vi=zeros(size(s,1),m,3);
  for i=1:2
    Vi(:,i,1)=funeval(cv{i},fspace{i},s);
    Vi(:,i,2)=funeval(cv{i},fspace{i},s,1);
    Vi(:,i,3)=funeval(cv{i},fspace{i},s,2);
  end
  V  =[Vi(s<sstar,1,1);Vi(s>=sstar,2,1)];
  Vs =[Vi(s<sstar,1,2);Vi(s>=sstar,2,2)];
  Vss=[Vi(s<sstar,1,3);Vi(s>=sstar,2,3)];


  close all
  figure(1)
  plot(s,V,'k',sstar,funeval(cv{1},fspace{1},sstar),'k*');
  title('Value Function');
  xlabel('S');
  ylabel('V');


  figure(2)
  plot(s,Vs,'k',sstar,funeval(cv{1},fspace{1},sstar,1),'k*');
  ylim([0 5])
  title('Marginal Value Function');
  xlabel('S');
  ylabel('V''');


  figure(3)
  plot(s,Vss,'k',sstar,funeval(cv{1},fspace{1},sstar,2),'k*');
  title('Curvature of Value Function');
  xlabel('S');
  ylabel('V"');
  ylim([-10 5])


  e=rho*V-H*(P-c./s).*s.*(s>sstar)-(alpha*(1-s)-(s>sstar)*H).*s.*Vs-0.5*sigma^2*s.^2.*Vss;
  figure(4)
  plot(s,e);
  title('Approximation Residual');
  xlabel('S');
  ylabel('Residual');


  disp(' ')
  disp('      S* ')
  disp([sstar])

prtfigs(mfilename,'Solution to the Fish Harvesting Model',[1 2 3 4])