% FUNDEFN Defines a function family structure
% USAGE
%   fspace = fundefn(bastype,n,a,b,order,s1,s2,...);
% INPUTS
%   bastype  : the string referencing function family ('cheb','spli' or 'lin')
%   n        : order of approximation along each dimension
%   a        : left endpoints of interpolation intervals in each dimension
%   b        : right endpoints of interpolation intervals in each dimension
%   order    : for 'spli' bastype, the order of the spline (default: 3 for cubic)
%   s1,s2... : additional column vectors for appending discrete variables to
%                the basis
%
% OUTPUT
%   fspace  : a function family structure
%
% USES: fundef
% See also: FUNDEF, FUNEVAL, FUNBAS.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function fspace = fundefn(bastype,n,a,b,order,varargin);

if ~isempty(bastype)
   d = length(n);
   if length(a) ~= d, error('a must be same dimension as n'); end
   if length(b) ~= d, error('a must be same dimension as n'); end
   if any(a>b), error('left endpoint must be less than right endpoint'); end
   if any(n<2), error('n(i) must be greater than one'); end
   if nargin<5 | isempty(order), order=3; end
   params = cell(1,d+length(varargin));
   switch bastype
      case 'cheb', for i=1:d; params(i)= {{'cheb',n(i),a(i),b(i)}}; end
      case 'spli', for i=1:d; params(i)= {{'spli',[a(i);b(i)],n(i)-order+1,order}}; end
      case 'lin',  for i=1:d; params(i)= {{'lin',[a(i);b(i)],n(i)}}; end   
      otherwise,   error ('basis type must be ''cheb'',''spli'', or ''lin''')
   end
else
   d = 0;
   params = cell(1,length(varargin));
end

for i=1:length(varargin)
  if all(abs(diff(diff(varargin{i})))<5e-15*mean(abs(varargin{i})))
    params(d+i) = {{'lin',varargin{i},length(varargin{i})}}; 
  else
    params(d+i) = {{'lin',varargin{i},0}}; 
  end
end
   
fspace = fundef(params{:});
