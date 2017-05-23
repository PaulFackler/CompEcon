% DEMDP07 Optimal Growth Model
  disp('DEMDP07 OPTIMAL GROWTH MODEL')
  close all  

% ENTER MODEL PARAMETERS
  alpha =  0.2;														% utility parameter
  beta  =  0.5;														% production elasticity
  gamma =  0.9;														% capital survival rate
  sigma =  0.1;														% production shock volatility
  delta =  0.9;														% discount factor

% COMPUTE SHOCK DISTRIBUTION
  m = 3;				      										% number of shocks
  [e,w] = qnwlogn(m,-sigma^2/2,sigma^2);           % shocks and proabilities

% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfdp07';						  		% model functions
  model.discount = delta;									% discount factor
  model.e = e;														% shocks
  model.w = w;														% probabilities
  model.params = {alpha beta gamma};	    % other parameters
  
% DEFINE APPROXIMATION SPACE
  n     = 10;	  													% degree of approximation
  smin  =  5;		  												% minimum state
  smax  = 10;			  											% maximum state
  fspace = fundefn('cheb',n,smin,smax);   % function space
  snodes = funnode(fspace);						    % state collocaton nodes
  
% COMPUTE CERTAINTY-EQUIVALENT STEADY-STATE
  estar = 1;				  				                         % certainty equivalent shock
  xstar = ((1-delta*gamma)/(delta*beta))^(1/(beta-1)); % steady-state action
  sstar = gamma*xstar + xstar^beta;							    	 % steady-state state
  pstar = (sstar-xstar).^(-alpha);                     % steady state shadow price  

% CHECK MODEL DERIVATIVES AT CE STEADY STATE
  dpcheck(model,sstar,xstar,estar);

% COMPUTE L-Q APPROXIMATION
  [vlq,xlq] = lqapprox(model,snodes,sstar,xstar,pstar);  
  
% SOLVE BELLMAN EQUATION
  [c,s,v,x,resid] = dpsolve(model,fspace,snodes,vlq,xlq);
   
% COMPUTE L-Q APPROXIMATION FOR PLOTTING 
  [vlq,xlq,plq] = lqapprox(model,s,sstar,xstar,pstar);
  
% PLOT OPTIMAL POLICY

  figure(1);
  plot(s,x./s,s,xlq./s,sstar,xstar/sstar,'*k');
  title('Optimal Investment Policy');
  legend('Chebyshev','L-Q');
  xlabel('Wealth');
  ylabel('Investment (% of Wealth)');

% PLOT VALUE FUNCTION
  figure(2);
  plot(s,v,s,vlq)
  title('Value Function');
  
  legend('Chebyshev','L-Q');
  xlabel('Wealth');
  ylabel('Value');

% PLOT SHADOW PRICE FUNCTION
  figure(3);
  
  p = funeval(c,fspace,s,1);
  plot(s,p,s,plq,sstar,pstar,'*k');
  title('Shadow Price Function');
  legend('Chebyshev','L-Q');
  xlabel('Wealth');
  ylabel('Price');

% PLOT RESIDUAL
  figure(4);
  plot(s,resid);
  title('Approximation Residual');
  xlabel('Wealth');
  ylabel('Residual');

% COMPUTE EXPECTED STATE AND POLICY PATH
  nyrs = 20; nrep = 2000;
  sinit = smin*ones(nrep,1);
  [spath,xpath] = dpsimul(model,sinit,nyrs,s,x);  

% PLOT EXPECTED STATE PATH
  figure(5);
  plot(0:nyrs,mean(spath));
  title('Expected Wealth');
  xlabel('Year');
  ylabel('Wealth');

% PLOT EXPECTED POLICY PATH
  figure(6);
  plot(0:nyrs,mean(xpath));
  title('Expected Investment');
  xlabel('Year');
  ylabel('Investment');
  
% COMPUTE STEADY-STATE DISTRIBUTION
  nsmooth = 5; nbin = 80;
  [ss,pi,xx] = dpstst(model,nsmooth,nbin,s,x);
  sstar = pi'*ss;  
  xstar = pi'*xx;

% PLOT STEADY-STATE DISTRIBUTION
  figure(7);
  h=bar(ss,pi,1); set(h,'FaceColor',[.75 .75 .75])
  title('Steady State Distribution');
  xlabel('Wealth');
  ylabel('Probability');
  fprintf('Steady State Means\n') 
  fprintf('   Wealth       = %5.2f\n'  ,sstar)
  fprintf('   Investment   = %5.2f\n'  ,xstar)
  
% SAVE PLOTS AS EPS FILES
  prtfigs(mfilename,'Solution to the Optimal Growth Model',[1 4 5 7])