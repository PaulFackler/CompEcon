%% LINBASE
%
%  Piecewise linear basis functions
%
%  Usage
%    [B,x] = linbase(breaks,evennum,x,order)
%  Input
%    breaks    : n.1 vector of breakpoints
%    a         : lower bound of approximation interval
%    b         : upper bound of approximation interval
%    x         : k.1 vector of the evaluation points (default: breaks)
%    order     : order of differentiation (default: 0)
%                if vector, cell array returned, otherwise matrix returned
%  Output
%    B         : k.n basis matrix or cell array of k.n basis matrices
%    x         : evaluation points (useful if defaults values are computed)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [B,x] = linbase(breaks,evennum,x,order)

if nargin<4 | isempty(order), order=0;
  if nargin<3, x = [];
    if nargin<2, evennum = 0;
      if nargin<1, error('At least one parameter must be passed');
      end; end; end; end;

% GET DEFAULTS
if isempty(x)
  x = linnode(breaks,evennum);
end
n = length(breaks);

% If multiple orders are requested make recursive call
% Inefficient but easy to code!
k = length(order);
if k>1
  B = cell(k,1);
  for ii=1:k
    B{ii} = linbase(breaks,evennum,x,order(ii));
  end
  return
end

if order~=0               % recursively generate differential operators
  [D,n,a,b,parms] = lindop(breaks,evennum,order);
  B = linbase(parms{:},x)*D{end};
  return
end

m = size(x,1);

% Determine the maximum index of
%   the breakpoints that are less than or equal to x,
%   (if x=b use the index of the next to last breakpoint).
if evennum
  ind = fix((x-breaks(1)).*((n-1)./(breaks(end)-breaks(1))))+1;
  ind = min(max(ind,1),n-1);
else
  ind = lookup(breaks,x,3);
end

z = (x-breaks(ind))./(breaks(ind+1)-breaks(ind));
B = sparse([1:m 1:m],[ind ind+1],[(1-z) z],m,n);


%% LINDOP
%
%  Differential operator matrices for a piecewise linear basis.
%
%  Usage
%    [D,n,a,b,parms] = lindop(breaks,evennum,order)
%  Input
%    breaks    : n.1 vector of breakpoints
%    evennum   : n if breakpoints are even, 0 otherwise
%    order     : desired order of differentiation (a scalar)
%  Output
%    D         : cell array of matrices of size abs(order) by 1.
%    n,a,b,parms : characteristics of the altered family of functions
%  Note: A piecewise linear function with n-1 pieces can be described by
%    the function values at each of the n breakpoints (f).  The derivative
%    is taken to be a piecewise linear function with values at the
%    breakpoints equal to the 3-point finite difference approximations at
%    the breakpoints. These values are produced by D{end}*f.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [D,n,a,b,parms] = lindop(breaks,evennum,order)

newbreaks = breaks;
n = length(breaks);
D = cell(abs(order),1);
for i=1:order
  d = 1./diff(newbreaks);
  d = sparse([1:n-1 1:n-1],[1:n-1 2:n],[-d d],n-1,n);
  if i>1
    D{i} = d*D{i-1};
  else
    D{1} = d;
  end
  newbreaks = (newbreaks(1:end-1)+newbreaks(2:end))/2;
  n = n-1;
end
for i=-1:-1:order
  newbreaks = [[3 -1]*newbreaks(1:2);
    (newbreaks(1:end-1)+newbreaks(2:end));
    [-1 3]*newbreaks(end-1:end)]/2;
  d = diff(newbreaks)';
  n = n+1;
  d = tril(d(ones(n,1),:),-1);
  if i<-1
    D{-i} = d*D{-i-1};
  else
    D{1} = d;
  end
  % adjustment to make value at original lower bound equal 0
  if evennum>0
    temp = linbase(newbreaks,length(newbreaks),breaks(1),0)*D{-i};
  else
    temp = linbase(newbreaks,0,breaks(1),0)*D{-i};
  end
  D{-i} = D{-i}-temp(ones(length(newbreaks),1),:);
end

if nargout>1
  n = length(newbreaks);
  a = newbreaks(1);
  b = newbreaks(end);
  if evennum>0
    parms = {newbreaks,n};
  else
    parms = {newbreaks,0};
  end
end