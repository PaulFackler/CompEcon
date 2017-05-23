% FUNEVAL2 Evaluates a function using 'direct'format
% USAGE
%   f=funeval2(g,B,order)
% INPUTS
%   B     : a basis structure in 'direct' format
%   g     : a coefficient matrix 
%   order : order of differentiation

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function f=funeval2(g,B,order)

  if nargin<3 | isempty(order)        % evaluate the function only
    kk=1;
    d=size(B.order,2);
    order=zeros(1,d);
  else
    [kk,d]=size(order);
  end

  % reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
  order=fliplr(order+ones(size(order,1),1)*(size(B.vals,1)*(0:d-1)-B.order+1));
  f=zeros(size(B.vals{1},1),size(g,2),kk);   % preallocate
  for i=1:kk
    f(:,:,i)=cdprodx(B.vals,g,order(i,:));
  end