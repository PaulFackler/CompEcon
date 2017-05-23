% DEMRE02 Asset Pricing Model
  
% ENTER MODEL PARAMETERS
  beta  = 0.5;                                    % coefficient of risk aversion
  mu    = [3.0 3.0 3.0 3];                              % long-run mean dividend
  gamma = [0.5 0.25 0.5 0.4];                              % dividend autoregression coefficient
  sigma = [0.05 0.1 0.1 0.05];                            % dividend volatility

  mu    = [3.0 3.0];                              % long-run mean dividend
  gamma = [0.5 0.25];                              % dividend autoregression coefficient
  sigma = [0.05 0.1];                            % dividend volatility

  rho=0.75;
  nassets=length(mu);
  C=rho+zeros(nassets,nassets);
  C((1:nassets)+(0:nassets:nassets*(nassets-1)))=1;
  delta = 0.9;                                    % discount factor

% COMPUTE SHOCK DISTRIBUTION
  m = 5;                                                % number of production shocks
  Sigma=diag(sigma)*C*diag(sigma);                      % shock covariance matrix  
  [e,w] = qnwnorm(m+zeros(1,nassets),zeros(1,nassets),Sigma); % normal nodes and weights

% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfre02';                                % model functions
  model.params = {delta mu gamma beta};                 % other parameters
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.explicit = 1;
  model.noxnext  = 0;
  
% DEFINE APPROXIMATION SPACE
  n      = 3;                                          % degree of approximation
  dmin   = mu+min(e)./(1-gamma);                       % minimum production
  dmax   = mu+max(e)./(1-gamma);                       % maximum production
  fspace = fundefn('cheb',n+zeros(1,nassets),dmin,dmax);    % function space
  dnodes  = funnode(fspace);                           % state collocaton nodes
  d=gridmake(dnodes);

    
% SOLVE RATIONAL EXPECTATIONS EQULIBRIUM
  nres=4;
  options=struct('usebroyden',1,...
                 'expapprox',1,...
                 'lowmemory',0,...
                 'tol', sqrt(eps),...
                 'stepsize',1.75,...
                 'maxit',500,...
                 'showiters',0,...
                 'nres',nres);

  % Solve the linearized problem
  N=size(d,1);
  dstar=mu;
  pstar=(delta/(1-delta))*mu;
  ddim=length(dstar);
  B22=sum(mu).^-beta*eye(ddim);
  B21= -beta/(1-delta)*sum(mu).^(-beta-1)*(mu'*ones(1,ddim));
  A22=delta*B22;
  A21=A22 - beta*delta*(2-delta)/(1-delta)*sum(mu).^(-beta-1)*(mu'*ones(1,ddim));
  P=diag(gamma);
  temp=B21-A21*P;
  C=reshape((kron(P,A22)-kron(eye(size(d,2)),B22))\temp(:),ddim,ddim);
  xinit=pstar(ones(size(d,1),1),:)+(d-dstar(ones(size(d,1),1),:))*C';
  u = sum(d,2).^(-beta);
  zinit  = u(:,ones(1,size(d,2))).*(xinit+d);
  if options.expapprox
    c=funfitxy(fspace,dnodes,zinit);
  else
    c=funfitxy(fspace,dnodes,xinit);
  end
tic
  [c,d,p,x,f] = resolve(model,fspace,options,c); 
toc
  
% PLOT EQUILIBRIUM PRICE FUNCTIONS  
  n=ones(1,length(d));
  for i=1:length(n)
    n(i)=length(d{i});
  end
  p=reshape(p,[n length(d)]);
  x=reshape(x,[n length(d)]);
  f=reshape(f,[n length(d)]);

  xx=cell(1,nassets);
  for i=1:nassets
    xx{i}=mu(i);
  end
  xx{1}=d{1};

  if options.expapprox
    z1=funeval(c,fspace,xx);
    p1=feval(model.func,'x',gridmake(xx),[],z1,[],[],[],model.params{:});
  else
    p1=funeval(c,fspace,xx);
  end

  figure(1);
  plot(d{1},p1);
  xlabel('d_1');
  legstr=repmat('p_1(d_1)',nassets,1);
  legstr(:,3)=num2str((1:nassets)');
  legend(legstr ,0)
  title('Equilibrium Price Functions (d_i=\mu_i for i>1)')
