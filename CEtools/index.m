% INDEX Converts between single and multiple indices
%         Replacement for SUB2IND and IND2SUB 
% USAGE
%   j=index(siz,ind)
% INPUTS
%   siz - 1xd vector of dimensions
%   ind - mxd or mx1 matrix of indices
% OUTPUT
%   j - mx1 matrix of single indices if ind is mxd 
%         (conversion from multiple to single) 
%     - mxd matrix of multiple indices if ind is mx1
%         (conversion from single to multiple) 
%
% This function is self-inverting, i.e., index(siz,index(siz,ind)) equals ind.
%
% Example:
% Consider a 3-dimensional array with 4, 5 and 3 elements associated with
% dimensions 1, 2 and 3. The relationship between the multiple and single 
% index forms can be represented as
%         
%      1  5  9 13 17     21 25 29 33 37    41 45 49 53 57
%      2  6 10 14 18     22 26 30 34 38    42 46 50 54 58 
%      3  7 11 15 19     23 27 31 35 39    43 47 51 55 59 
%      4  8 12 16 20     24 28 32 36 40    44 48 52 56 60
%
% where rows represent dimension 1, columns dimension 2 and the three 
% matrices represent dimension 3. 
% To determine the vectorized value of element (1,3,2) of the array use:
%    siz=[4 5 3]; ind=[1 3 2]; index(siz,ind)
% which returns 29. 
% To determine the multiple index of the 56th element in the array use:
%   siz=[4 5 3]; ind=56; index(siz,ind)
% which returns [4 4 3]. 

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function j=index(siz,ind)

  if nargin~=2
    error('Must pass two arguments');
  end
  [m,d]=size(siz);
  if m~=1
    error('SIZ must be a row vector');
  end

  [m,di]=size(ind);
 
  if d<di                                 % pads SIZ with ones
    siz=[siz ones(1,di-d)];
    d=di;
  end

  if d==1 & di==1                         % trivial case - both 1-D
    if max(ind,[],1)>siz
      error('Index out of bounds');
    end
    j=ind;
  elseif di==d                            % convert multiple to single
    k=cumprod([1;siz(:)]);
    if any(max(ind,[],1)>siz)
      error('Index out of bounds');
    end
    j=(ind-1)*k(1:d)+1;
  elseif di==1                            % convert single to multiple
    k=cumprod([1;siz(:)]);
    if any(ind>k(d+1))
      error('Index out of bounds');
    end
    ind=ind-1;
    j=zeros(m,d);
    for i=d:-1:2
      j(:,i)=fix(ind./k(i));
      ind=rem(ind,k(i));
    end;
    j(:,1)=ind;
    j=j+1;
  else
    error('Arguments are column incompatible');
  end