% FOURBAS Defines basis matrices for Fourier series
% USAGE
%   x=fourbas(m,rho,beta,x,order);
% INPUTS
%   m       : the number of cos and sin terms to include
%             (there will be 2m+1 basis functions)
%   rho     : the periodicity of the function (default: 1)
%   beta    : exponent on the first term (default: 0)
%   x       : k-vector of the evaluation points (default: evenly spaced)
%   order   : the order of differentiation (default: 0)
%             if a vector, SPLIBAS returns a cell array 
%             otherwise it returns a matrix
% OUTPUTS
%   B :  a kxn basis matrix or cell array of basis matrices
%   x :  evaluation points (useful if defaults values are computed)
%
% Useful in fitting functions that are known to be periodic
%   with period rho, i.e. f(z+rho)=f(z).
%
% See also: FOURNODE, FOURDOP, FUNBAS, FUNEVAL.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [B,x]=fourbas(m,rho,beta,x,order);

  if nargin<3, error('Three parameters must be specified'), end
  if nargin<4, x=[]; end
  if nargin<5 | isempty(order), order=0; end
  
  if isempty(x)
    x=fournode(m,rho,beta);
  end 

  % If multiple orders are requested make recursive call
  % Inefficient but easy to code!
  k=length(order);
  if k>1
    B=cell(k,1);
    for ii=1:k
      B{ii}=fourbas(m,rho,beta,x,order(ii));
    end
    return
  end
  
  if order~=0               % recursively generate differential operators
    [D,n,parms]=fourdop(m,rho,beta,order);
    B=fourbas(parms{:},x)*D;     
    return
  end    

  n=m+m+1;
  B=zeros(size(x,1),n);
  nc=fix((n+1)/2);
  ns=fix(n/2);
  B(:,1)=x.^beta;
  z=(2*pi/rho)*x;       % normalize to [0,2pi]
  if n>1
    B(:,[2:2:2*ns])=sin(z*(1:ns));
  end
  if n>0
    B(:,[1:2:2*nc-1])=cos(z*(0:nc-1));
  end

    










