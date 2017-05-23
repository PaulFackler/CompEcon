% FISHH Solution function for fish harvesting problem
% See DEMFB03 for demonstration
function [sstar,cl,cdefl,cu,cdefu]=fishh(rho,alpha,sigma,H)

% Set up basis matrices
  a=log(0.005);        % lower bound
  b=log(10);           % upper bound
  nl=25;   nu=15;      % number of nodes for functions
  nbl=2;   nbu=1;      % number of boundary constraints on functions
  cdefl=fundef({'cheb',nl-nbl,a,0});
  cdefu=fundef({'cheb',nu-nbu,0,b});
  yl=funnode(cdefl);   
  yu=funnode(cdefu);   
  cdefl=fundef({'cheb',nl,a,0});
  cdefu=fundef({'cheb',nu,0,b});
  eyl=exp(yl);         
  eyu=exp(yu);
  Dl=funbas(cdefl,yl,1);   
  Du=funbas(cdefu,yu,1);
  B=rho*funbas(cdefl,yl)...
       -(alpha-0.5*sigma.^2)*Dl...
       -(0.5*sigma.^2)*funbas(cdefl,yl,2);
  temp=rho*funbas(cdefu,yu)...
       -(alpha-0.5*sigma.^2-H)*Du...
       -(0.5*sigma.^2)*funbas(cdefu,yu,2);
  B=[B zeros(nl-nbl,nu);zeros(nu-nbu,nl) temp];
% Add boundary constraints
  B=[B; ...
     funbas(cdefl,0)   -funbas(cdefu,0); ...          % V continuous at y=0
     funbas(cdefl,0,1) -funbas(cdefu,0,1)];           % Vx continuous at y=0
  if nbl==2; B=[B;funbas(cdefl,a,2) zeros(1,nu)]; end % lower boundary 
  if nbu==2; B=[B;zeros(1,nl) funbas(cdefu,b,2)]; end % upper boundary
% Basis for Vy
  D=[alpha*eyl*ones(1,nl).*Dl    zeros(nl-nbl,nu); ...
     zeros(nu-nbu,nl)            alpha*eyu*ones(1,nu).*Du; ...
     zeros(nbl+nbu,nl+nu)];
% RHS of DE residual function
  f=[zeros(nl-nbl,1);H*eyu;zeros(nbl+nbu,1)];
% Basis for residual function (Vy(0)=S*)
  phil10=funbas(cdefl,0,1);
% find the cutoff stock level  
  sstar=broyden('fishhr',0.5*(1-rho/alpha),[],B,D,f,phil10,nl);
% Break apart the coefficient vector and create structures to return
  c=(B+sstar*D)\(sstar*f);
  cl=c(1:nl);
  cu=c(nl+1:end);