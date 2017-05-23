% FUNEVAL1 Evaluates a function using 'tensor'format
% USAGE
%   f=funeval1(g,B,order);
% INPUTS
%   g     : a coefficient matrix 
%   B     : a basis structure in 'tensor' format
%   order : order of differentiation

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function f=funeval1(g,B,order)

  if nargin<3 | isempty(order)        % evaluate the function only
    kk=1;
    d=size(B.order,2);
    order=zeros(1,d);
  else                                % evaluate according to order
    [kk,d]=size(order);
  end
  
  % reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
  order=fliplr(order+ones(size(order,1),1)*(size(B.vals,1)*(0:d-1)-B.order+1));
  m=1; for j=1:d, m=m*size(B.vals{1,j},1); end
  f=zeros(m,size(g,2),kk);          % preallocate output matrix
  for i=1:kk
    f(:,:,i)=ckronx(B.vals,g,order(i,:));
  end