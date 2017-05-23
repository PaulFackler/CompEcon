% DEMSC02 Renewable Resource Management Model
function demsc02
  close all
  optset('scsolve','defaults');

  disp(' ')
  disp('DEMSC02 Renewable Resource Management Model')

% Set parameter values
  alpha = 0.5;
  beta  = 1;    
  K     = 1;   
  b     = 1;
  eta   = 0.5;
  c     = 5;
  gamma = 2;
  sigma = 0.1;
  rho   = 0.05;  

% 3 cases with known solutions 
  casenum=1;
  switch casenum
  case 1
    beta=1;    eta=0.5;  gamma=2;
  case 2
    beta=0;    eta=1;    gamma=1;
  case 3
    beta=-0.5; eta=2;    gamma=0.5;  % should increase smax with this option
  end

  clear model
  model.func='mfsc02';
  model.params={alpha,beta,K,b,eta,c,gamma,sigma,rho};

% Define the approximating family and nodes
  n = 35;
  smin = 0.2;
  smax = 1.2;
  fspace = fundefn('cheb',n,smin,smax);
  s = funnode(fspace);

% Call the stochastic control solver
  [cv,S,v,x,resid] = scsolve(model,fspace,s);

% Compute the known solution for this model at many values of S
  switch casenum
  case 1
    phi=2*(b/(alpha+rho-sigma.^2)).^2* ...
       (1+sqrt(1+c*((alpha+rho-sigma.^2)/b).^2));
    V=-phi*(1./S+alpha/(rho*K));
    Vs=phi./S.^2;
    q=b/sqrt(c+phi)*S;
  case 2
    temp=(alpha+rho)/(b+c*(alpha+rho));
    phi=log(b*temp)-c*temp+(alpha*log(K)-b*temp-0.5*sigma.^2)/(alpha+rho);
    phi=phi*b/rho;
    V=b/(alpha+rho)*log(S)+phi;
    Vs=b./((alpha+rho)*S);
    q=b*temp*S;
  case 3
    phi=c-sqrt(c.^2+2*b/(alpha+rho+sigma.^2/8));
    V=-phi*(sqrt(S)+alpha*sqrt(K)/rho);
    Vs=-0.5*phi./sqrt(S);
    q=b/(c-0.5*phi).^2*S; 
  end

% Plot the errors in the marginal value function and optimal control
  Vsa=funeval(cv,fspace,S,1);
  qa=feval(model.func,'x',S,[],Vsa,model.params{:});
  ee=[Vs-Vsa 100*(q-qa)];

  figure(1);  
  plot(S,ee(:,1),'-',S,ee(:,2),'--')
  title('Approximation Errors')
  xlabel('S')
  ylabel('Error')
  legend('Marginal Value Function','Optimal Control (x100)',2)

  figure(2);  
  plot(S,funeval(cv,fspace,S))
  title('Value Function')
  xlabel('S')
  ylabel('V')

  figure(3);  
  plot(S,Vsa)
  title('Shadow Price')
  xlabel('S')
  ylabel('V''(S)')

  figure(4);  
  plot(S,qa)
  title('Optimal Control')
  xlabel('S')
  ylabel('q')

  figure(5);
  plot(S,resid);
  title('Approximation Residual')
  xlabel('S')
  ylabel('Residual')

  [cp,Ex]=itodensity(model,fspace,cv);
  figure(6)
  plot(S,funeval(cp,fspace,S))
  title('Long-Run State Density')
  xlabel('S')
  ylabel('Probability')
  % set lower y limit to 0
  yy=get(gca,'ylim');yy(1)=0; set(gca,'ylim',yy); 

prtfigs(mfilename,'Solution to the Renewable Resource Model',[3 6 5 1])

  