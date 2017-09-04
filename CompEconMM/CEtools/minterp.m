%% MINTERP
%
%  Multidimensional linear interpolation
%
%  Usage
%    y = minterp(X,Y,x,even)
%  Input
%    X         : 1.d cell array of coordinate vectors
%                each must be sorted in ascending order
%                may be vector rather than a cell array for d=1
%    Y         : function values at grid points defined by X (use 
%                gridmake to expand gridpoints)
%                n.p array where n is the total number of grid
%                points defined by X (# of rows in gridmake(X))
%    x         : values at which to interpolate, x can be either:
%                  1.d cell array of vectors to produce interpolates on a new grid
%                  k.d matrix to produce interpolates at k arbitrary points
%    even      : 1 if vectors in X evenly spaced (speeds up algorithm if true)
%  Output
%    y         : interpolated values
%                N.p if x is a cell array defining a grid with N points
%                k.p is x is a kxd matrix

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function y = minterp(X,Y,x,even)

if nargin<4 || isempty(even)
  even = 0;
end

d = size(x,2);
if d==1
  if isa(X,'cell')
    X = X{1};
  end
  if even
    ind = evenlookup(X,x);
  else
    ind = lookup(X,x,3);
  end
  z = X(ind);
  y = Y(ind,:);
  ind = ind+1;
  z = (x-z)./(X(ind)-z);
  m = size(Y,2);
  if m>1, z = z(:,ones(size(Y,2),1)); end
  y = y.*(1-z)+Y(ind,:).*z;
else
  if isa(x,'cell')       % use tensor products
    B = cell(1,d);
    for i=1:d
      B{i} = getbas(X{i},x{i},even);
    end
    y = ckronx(B,Y,d:-1:1);
  else                    % use direct (row-wise tensor) products
    B = cell(1,d);
    for i=1:d
      B{i} = getbas(X{i},x(:,i),even);
    end
    y = cdprodx(B,Y,d:-1:1);
  end
end

function B = getbas(X,x,even)
m = size(x,1);
n = size(X,1);
if even
  ind = evenlookup(X,x);
else
  ind = lookup(X,x,3);
end
z = (x-X(ind))./(X(ind+1)-X(ind));
B = sparse([1:m 1:m],[ind ind+1],[(1-z) z],m,n);

function ind = evenlookup(X,x)
n = length(X);
ind = ceil((x-X(1))*(n/(X(end)-X(1))));
ind = min(max(ind,1),n-1);