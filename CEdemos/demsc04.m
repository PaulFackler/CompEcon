% DEMSC04 Optimal Fish Harvest Model
  close all  
  optset('scsolve','defaults');
  
  disp(' ')
  disp('DEMSC04 Optimal Fish Harvest Model')

% Define problem parameters
  alpha = 0.5;  
  sigma = 0.5;
  H     = 1; 
  P     = 1; 
  c     = 0.25;  
  rho   = 0.1;

% Pack model structure
  clear model
  model.func='mfsc04';
  model.params={alpha,sigma,H,P,c,rho};
 
% Define nodes and basis matrices
  n=150;  
  smin=0.01;
  smax=1.5;
  fspace=fundefn('lin',n,smin,smax);
  s=funnode(fspace);

% Define initial values
  v0=sqrt(s);
  v0=10*log(s+1);
  
% Call solver
  [cv,s,V,x,resid]=scsolve(model,fspace,s,v0);

% approximate the optimal switch point and value functions
  S=linspace(smin,1,501)';
  Vs=funeval(cv,fspace,S,1);
  x=feval(model.func,'x',S,[],Vs,model.params{:});
  i=sum(x==0);
  sstar=(S(i)+S(i+1))/2;
  S=[S(1:i);sstar;S(i+1:end)];
  Vs=[Vs(1:i);funeval(cv,fspace,sstar,1);Vs(i+1:end)];
  V=funeval(cv,fspace,S);
  Vss=funeval(cv,fspace,S,2);
 
% Produce plots
  close all
  figure(1)
  plot(S,V,'k',sstar,funeval(cv,fspace,sstar),'k*');
  title('Value Function');
  xlabel('S');
  ylabel('V');

  figure(2)
  plot(S,Vs,'k',sstar,funeval(cv,fspace,sstar,1),'k*',S,P-c./S,'r--');
  ylim([0 5])
  title('Marginal Value Function');
  xlabel('S');
  ylabel('V''');

  figure(3)
  plot(S,Vss,'k',sstar,funeval(cv,fspace,sstar,2),'k*');
  title('Curvature of the Value Function');
  xlabel('S');
  ylabel('V"');
  ylim([-10 5])

  [cp,Ex]=itodensity(model,fspace,cv);
  p=funeval(cp,fspace,s);
  fspace0=fundefn('cheb',51,fspace.a,3);
  s0=linspace(fspace0.a,fspace0.b,201)';
  cv0=funfitxy(fspace0,s0,(20*P)*s0);
  [cp0,Ex0]=itodensity(model,fspace0,cv0);
  p0=funeval(cp0,fspace0,s0);
  figure(4)
  plot(s,p,s0,p0)
  title('Long-run Density of Fish Stocks');
  xlabel('S')
  ylabel('Probability')
  % set lower y limit to 0
  yy=get(gca,'ylim');yy(1)=0; set(gca,'ylim',yy); 
  legend({'With harvesting','Without harvesting'})

  figure(5)
  plot(s,resid)
  title('Approximation Residual')
  xlabel('S')
  ylabel('Residual')

  disp(' ')
  disp('      S*       E[S]    E[S|x=0]')
  disp([sstar Ex Ex0])
  disp(' ')
  disp('Percentage of time inactive')
  disp(funeval(cp,fspace,sstar,-1))

prtfigs(mfilename,'Solution to the Fish Harvesting Model',[1 2 3 4])