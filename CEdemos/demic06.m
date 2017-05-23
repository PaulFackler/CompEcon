% DEMIC06 Optimal Fish Harvest Model (infinite harvest rate)
  close all
  disp(' ')
  disp('DEMIC06 Optimal Fish Harvest Model')

% Define parameters
  alpha = 0.5;  
  sigma = 0.5;
  P     = 1; 
  c     = 0.25;  
  rho   = 0.1;

% Create model variable
  clear model
  model.func='mfic06';
  model.params={alpha,sigma,P,c,rho};
  model.xindex=[0 0;2 1];
  model.F=[0;0];

% Define starting values 
  x=[0.01 0;1 0];

% Call solver
  n=50;  
  [cv,fspace,x]=icsolve(model,x,n,'cheb');
  sstar=x(2,1);

% Plot results
  s=linspace(min(0.01,x(1,1)),1.5,101)';
  s=sort([s;sstar;sstar+sqrt(eps)]);

  V=funeval(cv,fspace,s);
  Vs=funeval(cv,fspace,s,1);
  Vss=funeval(cv,fspace,s,2);
 
  ind=s>sstar;
  V(ind)=P*(s(ind)-sstar)-c*log(s(ind)/sstar)+funeval(cv,fspace,sstar);
  Vs(ind)=P-c./s(ind);
  Vss(ind)=c./s(ind).^2;

  figure(1)
  plot(s,V,'k',sstar,funeval(cv,fspace,sstar),'k*');
  title('Value Function');
  xlabel('S');
  ylabel('V');
  xlim([s(1) s(end)])

  figure(2)
  plot(s,P-c./s,'r:',s,Vs,'k',sstar,funeval(cv,fspace,sstar,1),'k*');
  ylim([0 5])
  title('Marginal Value Function');
  xlabel('S');
  ylabel('V''');
  xlim([s(1) s(end)])

  figure(3)
  plot(s,Vss,'k',sstar,funeval(cv,fspace,sstar,2),'k*');
  title('Curvature of Value Function');
  xlabel('S');
  ylabel('V"');  
  xlim([s(1) s(end)])
  ylim([-10 5])

  % Long run density of fish stocks
  [cp,Ex]=itodensity(model,fspace);
  p=funeval(cp,fspace,s);
  p(s>sstar)=0;
  fspace0=fundefn('cheb',51,fspace.a,3);
  s0=linspace(fspace0.a,fspace0.b,201)';
  cv0=funfitxy(fspace0,s0,(20*P)*s0);
  [cp0,Ex0]=itodensity(model,fspace0);
  p0=funeval(cp0,fspace0,s0);
  figure(4)
  plot(s,p,s0,p0)
  title('Long-run Density of Fish Stocks');
  xlabel('S')
  ylabel('Probability')
  % set lower y limit to 0
  yy=get(gca,'ylim');yy(1)=0; set(gca,'ylim',yy); 
  legend({'With harvesting','Without harvesting'})

  disp('      S*       E[S]')
  disp([sstar Ex])

  prtfigs(mfilename,'Solution to the Fish Harvesting Model',[1 2 3 4])