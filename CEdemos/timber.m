% Timber An alternative way to solve the timber managment problem
%   [sstar,s,v,p]=timber([],model);
% See Also: demdp02
function [out1,out2,out3,out4]=timber(sstar,model,z,Phi,phi0,phi1,price,C,cc,alpha,sigma)
if nargin>2                    % compute residuals
  delta=model.discount;
  phi0=delta*ones(size(z))*phi0;
  B=zeros(size(z,1));
  b=zeros(size(z));
  snext=z*sstar;
  ind=ones(size(z));
  dk=delta;
  for i=1:50
    if all(snext>sstar), break; end
    snext=cc*(1-exp(-alpha))+exp(-alpha)*snext;
    ii = snext>sstar & ind==1;
    B(ii,:)=dk*phi0(ii,:);
    b(ii)=dk*(price*snext(ii)-C);
    dk=delta*dk;
    ind=ind-ii;
  end
  B=Phi-B;
  c=B\b;
  res=price*sstar-C+(phi0(1,:)-phi1)*c;
  out1=res;
  out2=c;
else                            % set up collocation and call solver
  % Define function space, nodes and basis matrices
  n     =  100;                           % degree of approximation
  fspace = fundefn('lin',n,0,1);          % function space
  z = funnode(fspace);                    % state collocaton nodes
  Phi=funbas(fspace,z);
  phi0=funbas(fspace,0);
  phi1=funbas(fspace,1);

  % call root finder
  cc=model.params{3};            % use as initial condition
  sstar=broyden(mfilename,cc/2,[],model,z,Phi,phi0,phi1,model.params{:});
  [res,c]=feval(mfilename,sstar,model,z,Phi,phi0,phi1,model.params{:});

  % compute solute values
  z=linspace(0,1,501)';
  s=z*sstar;
  v=funeval(c,fspace,z);
  p=funeval(c,fspace,z,1)/sstar;

  out1=sstar; out2=s; out3=v; out4=p;
end