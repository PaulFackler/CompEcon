% FUNEVAL3 Evaluates a function using 'expanded'format
% USAGE 
%   f=funeval3(g,B,order);
% INPUTS
%   B     : a basis structure in 'expanded' format
%   g     : a coefficient matrix 
%   order : order of differentiation

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function f=funeval3(g,B,order)

  if nargin<3 | isempty(order)   % default: evaluate the function only
    if isa(B.vals,'cell')
      kk=length(B.vals);
      order=1:kk;
    else
      kk=1;
      order=1;
    end
  else                           % determine which bases are to be used 
    kk=size(order,1);
  end

  if isa(B.vals,'cell')
    m=size(B.vals{1},1);
    f=zeros(m,size(g,2),kk); 
    for i=1:kk
      % Determine which element of B.vals is the desired basis, if any
      ii=find(ismember(B.order,order(i,:),'rows'));
      if isempty(ii)
        error('Requested basis matrix is not available');
      end
      if length(ii)>1,
        warning('Redundant request in FUNEVAL3')
        ii=ii(1);
      end
      f(:,:,i)=B.vals{ii}*g;
    end
  else
    m=size(B.vals,1);
    f=zeros(m,size(g,2),kk); 
    for i=1:kk
      f(:,:,i)=B.vals*g;
    end
  end 