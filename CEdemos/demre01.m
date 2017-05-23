% DEMRE01 Optimal Growth Model
  close all  

% ENTER MODEL PARAMETERS
  alpha = 3;                       % utility parameter
  beta  = 0.33;                    % production elasticity
  gamma = 1;                       % capital survival rate
  rho   = 0.9;                     % production shock autocorrelation
  delta = 0.95;                    % discount factor
  sigma = 0.02;                    % production shock volatility

% COMPUTE SHOCK DISTRIBUTION
  K = 3;                                           % number of shocks
  [e,w] = qnwnorm(K,0,sigma^2);                    % shocks and proabilities

% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfre01';                          % model function file
  model.e = e;                                     % shocks
  model.w = w;                                     % probabilities
  model.params = {alpha beta gamma rho delta};     % other parameters
  model.explicit = 1;

% COMPUTE CERTAINTY-EQUIVALENT STEADY-STATE
  sstar = [((1/delta-gamma)/beta)^(1/(beta-1)) 0];
  xstar = sstar(1)^beta+(gamma-1)*sstar(1);
  zstar = feval(model.func,'h',sstar,xstar,[],w'*e,sstar,xstar,model.params{:});
  [ss,xx,zz]=ress(model,sstar,xstar);
% COMPUTE LINEARIZED SOLUTION
  [C,P,D,A,B]=relin(model,sstar,xstar);
  
% DEFINE APPROXIMATION SPACE
  n     = [10 6];                             % degree of approximation
  smin  = [sstar(1)/2    min(e)/(1-rho)];     % minimum state
  smax  = [1.5*sstar(1)  max(e)/(1-rho)];     % maximum state
  fspace = fundefn('cheb',n,smin,smax);       % function space
  snodes = funnode(fspace);                   % state collocaton nodes
  S=gridmake(snodes);
 
% Use linearized solution as starting values
  xinit=xstar+(S-sstar(ones(size(S,1),1),:))*C';
  
  snext=feval(model.func,'g',S,xinit,[],0,[],[],model.params{:});
  xnext=xstar+(snext-sstar(ones(size(S,1),1),:))*C';
  zinit=feval(model.func,'h',S,xinit,[],0,snext,xnext,model.params{:});

% Set algorithm options  
  options=struct('usebroyden',  1,...
                 'expapprox',   0,...
                 'lowmemory',   0,...
                 'maxit',     500,...
                 'stepsize', 1.75,...
                 'showiters',   1,...
                 'checks',      0,...
                 'nres',        5);

  if options.expapprox
    c=funfitxy(fspace,snodes,zinit);
  else
    c=funfitxy(fspace,snodes,xinit);
  end

  [c,s,x,z,f,resid] = resolve(model,fspace,options,c);


  K=s{1};
  V=s{2};
  x=reshape(x,size(K,1),size(V,1));
  z=reshape(z,size(K,1),size(V,1));
  f=reshape(f,size(K,1),size(V,1));
  resid=reshape(resid,size(K,1),size(V,1));
  ind=[1 9 16 23 31];
  x=x(:,ind);
  z=z(:,ind);
  f=f(:,ind);
  resid=resid(:,ind);

  [KK VV]=gridmake(K-sstar(1),V(ind)-sstar(2));


% PLOT OPTIMAL POLICY
  figure(1);
  plot(K,x)
  title('Optimal Consumption Policy');
  xlabel('Wealth');
  ylabel('Consumption');
  legend([repmat('V=',5,1) num2str(V(ind))],0)

% PLOT OPTIMAL EXPECTATION FUNCTION
  figure(2);
  plot(K,z)
  title('Optimal Expectation Function');
  xlabel('Wealth');
  ylabel('z');
  legend([repmat('V=',5,1) num2str(V(ind))],0)
  
% PLOT EQUILIBRIUM FUNCTION
  figure(3);
  plot(K,f);
  title('Equilibrium Function');
  xlabel('Wealth');
  legend([repmat('V=',5,1) num2str(V(ind))],0)

% PLOT RESIDUAL
  figure(4);
  plot(K,resid);
  title('Approximation Residual');
  xlabel('Wealth');
  ylabel('Residual');
  legend([repmat('V=',5,1) num2str(V(ind))],0)
