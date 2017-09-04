%% FUNDEFN 
%
%   Defines a basis function family and computes associated standard
%   interpolation matrix and nodes.
%
% Usage
%   [basis,Phi,xnodes,xcoord] = fundefn(basistype,n,a,b,order)
% Let
%   d = dimension of the basis function domain
%   N = prod(n) number of basis functions
% Input
%   basistype : string indicating basis function type ('cheb','spli' or 'lin')
%   n         : 1.d ordder of approximation per dimension
%   a         : 1.d lower bounds of approximation interval per dimension
%   b         : 1.d upper bounds of approximation interval per dimension
%   order     : for 'spli' basistype, the order of the spline (default: 3 for cubic)
% Output
%   basis     : function family structure
%   Phi       : N.N interpolation matrix
%   xnodes    : N.d interpolation nodes
%   xcoord    : 1.d cell array of node coordinates by dimension

% Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [basis,Phi,xnodes,xcoord] = fundefn(basistype,n,a,b,order)

d = length(n);
if length(a) ~= d, error('a must be same dimension as n'); end
if length(b) ~= d, error('a must be same dimension as n'); end
if any(a>b), error('lower bound must be less than upper bound'); end
if any(n<2), error('n(i) must be greater than one'); end
if nargin<5 || isempty(order), order=3; end

% params = cell(1,d);
% switch basistype
%   case 'cheb', for i=1:d; params(i)= {{'cheb',n(i),a(i),b(i)}}; end
%   case 'spli', for i=1:d; params(i)= {{'spli',[a(i);b(i)],n(i)-order+1,order}}; end
%   case 'lin',  for i=1:d; params(i)= {{'lin',[a(i);b(i)],n(i)}}; end
%   otherwise,   error ('basis type must be ''cheb'',''spli'', or ''lin''')
% end
% basis = fundef(params{:});

params = cell(1,d);
if iscell(basistype)
  for i=1:d
    switch basistype{i}
      case 'cheb', params(i) = {{'cheb',n(i),a(i),b(i)}};
      case 'spli', params(i) = {{'spli',[a(i);b(i)],n(i)-order+1,order}};
      case 'lin',  params(i) = {{'lin',[a(i);b(i)],n(i)}};
      otherwise,   error ('basis type must be ''cheb'',''spli'', or ''lin''')
    end
  end
else
  for i=1:d
    switch basistype
      case 'cheb', params(i) = {{'cheb',n(i),a(i),b(i)}};
      case 'spli', params(i) = {{'spli',[a(i);b(i)],n(i)-order+1,order}};
      case 'lin',  params(i) = {{'lin',[a(i);b(i)],n(i)}};
      otherwise,   error ('basis type must be ''cheb'',''spli'', or ''lin''')
    end
  end
end
basis = fundef(params{:});

if nargout>1, Phi = funbase(basis); end
if nargout>2, [xnodes,xcoord] = funnode(basis); end



%% FUNDEF
%
%  Creates a fauction family definition structure
%
%  Usage
%      g = fundef({bastype1,p11,p12,...},{bastype2,p21,p22,...},...);
%    or
%      [g,errorstr] = fundef({bastype1,p11,p12,...},{bastype2,p21,p22,...},...);
%    First syntax will generate error message and halt execution if error
%    encountered. Second syntax sets g=[] and returns error message if an
%    error in encountered, otherwise it returns a function definition
%    structure and an empty string (test using if isempty(errorstr)).
%  Input
%    fundef accepts any number of input arguments, each of which must be a
%    cell array containing the following information:
%      bastype     : string prefix referencing function family name ('cheb' or 'spli')
%      p1, p2, ... : parameters required by the specified family
%    The parameters required by a specific family are documented in that
%    family's DEF file. For example the CHEB family requires n, the degree
%    of the polynomial (1 plus its order), and a and b, the bounds.
%  Examples:
%    cdef = funbase({'cheb',10,0,1});
%      defines a 1-D polynomial of order 9 on the interval [0,1].
%    cdef = funbase({'cheb',10,0,1},{'cheb',5,-1,1});
%      defines a 2-D polynomial whose first dimension is a polynomial of
%      order 9 on the interval [0,1] and whose second dimension is a
%      polynomial of order 4 on the interval [-1 1].
%    cdef = funbase({'cheb',10,0,1},{'spli',[-1;1],5});
%      defines a 2-D function whose first dimension is order 9 polynomial
%      on the interval [0,1] and whose second dimension is a cubic spline
%      with 5 evenly spaced breakpoints on the interval [-1 1].

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [g,errorstr] = fundef(varargin)

errorstr = '';
g = [];
% !!!!!!!!!! perform consistency checks !!!!!!!!!!!!!!
if isstruct(varargin{1})
  if nargin>1
    errorstr = 'Pass a single function structure for consistency check';
    if nargout<2, error(errorstr); else return; end
  end
  c = varargin{1};
  f = getfields(c);
  if ~strcmp(f,{'d';'n';'bastype';'parms'})
    errorstr = 'Structure has improper fields';
    if nargout<2, error(errorstr); else return; end
  end
  [m,d] = size(c.n);
  if m~=1 || d~=c.d
    errorstr = '''n'' field of improper size';
    if nargout<2, error(errorstr); else return; end
  end
  [m,d] = size(c.bastype);
  if m~=1 || d~=c.d
    errorstr = '''bastype'' field of improper size';
    if nargout<2, error(errorstr); else return; end
  end
  [m,d] = size(c.parms);
  if m~=1 || d~=c.d
    errorstr = 'parameter field of improper size';
    if nargout<2, error(errorstr); else return; end
  end
  if size(c.vals,1)~=prod(c.n)
    errorstr = 'Rows in ''vals'' field does not equal product of ''n'' field';
    if nargout<2, error(errorstr); else return; end
  end
  if ~iscell('c.bastype')
    errorstr = '''bastype'' field should be 1xd cell array (cell string)';
    if nargout<2, error(errorstr); else return; end
  end
  if ~iscell('c.parms')
    errorstr = '''parms'' field should be 1xd cell array';
    if nargout<2, error(errorstr); else return; end
  end
  g = c;  % no errors found
  
else  % !!!!!!!! Create a new coefficient structure !!!!!!!!
  
  errorstr = '';
  g = [];
  if isa(varargin{1},'cell')
    d = length(varargin);
  else d = 1;
    varargin(1) = varargin{:};
  end
  n = zeros(1,d);
  b = zeros(1,d);
  a = zeros(1,d);
  parms = cell(1,d);
  
  for j=1:d
    bastype{j} = varargin{j}{1};
    varargin{j}(1) = [];
    if ~exist([bastype{j} 'base'])
      errorstr = ['Cannot find M-file ' lower(bastype{j}) 'base.m'];
      if nargout<2
        error(errorstr)
      else
        return
      end
    end
  end
  
  warning = [];
  for j=1:d
    auxname = [bastype{j} 'def'];
    if exist(auxname,'file');  % check if DEF file exists
      m = nargin(lower(auxname));
      if m>0 && m<length(varargin{j})
        errorstr = ['Too many parameters for ' lower(bastype{j}) 'def'];
        if nargout<2
          error(errorstr);
        else
          return
        end
      end
      eval(['[n(j),a(j),b(j),parms{j}] = ' bastype{j} 'def(varargin{j}{:});']);%,'errorstr=lasterr;');
      if ~isempty(errorstr)
        if nargout<2
          error(errorstr);
        else
          return
        end
      end
    else
      warning = {warning;upper(auxname)};
    end
  end
  if ~isempty(warning)
    disp('Error:')
    disp('The following parameter definition files could not be found: ')
    disp(warning{2:end});
    error(' ');
  end
  
  % Assign values to fields
  g = struct('d',[],'n',[],'a',[],'b',[],'bastype',[],'parms',[]);
  g.d = d;
  g.n = n;
  g.a = a;
  g.b = b;
  g.bastype = bastype;
  g.parms = parms;
  
end