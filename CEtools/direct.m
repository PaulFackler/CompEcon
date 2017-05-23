% DIRECT DIvided RECTangles global optimization algorithm
% Solves
%    min_x f(x) s.t. a <= x <= b   
% USAGE
%   [xbest,fbest,convtype,z] = direct(f,bounds,options,varargin);
% INPUTS
%   f        : objective function (string or function handle)
%   bounds   : a dx2 matrix of bounds: bounds=[a b]
%   options  : a structure variable (described below)
%   varargin : additional parameters to pass to f
% OUTPUTS
%   xbest    : the optimal point
%   fbest    : the value at the optimal point
%   z        : a structure containing a complete record of the algorithm
%                (described below)
%
% Fields in options structure
%   maxflag    set to 1 for max problems, 0 for min problems [0]
%   sizeconst  constant used in rectangle size function  [0.5]
%   distopt    1/0 for distance/volume measure of size [1]
%   oneeach    1 if only 1 rectangle of each size is split per iteration [0]
%   greedy     if 1, use all sizes greater than the min slope size [0]
%   epsilon    global/local weight parameter [1e-4]
%   maxit      maximum number of iterations [50]
%   maxeval    maximum # of function evaluations [10000]
%   maxcuts    maximum # of side lengths [32]
%   fopt       optimal value of the objective function [-inf]
%   tol        relative error tolerance parameter [0.0001]
%   minlength  stop if best rectangle has all sides 1ess than this [1e-4]
%   minevals   but must evaluate at least this many points [0]
%   ShowIters  prints iteration results [1]
%
% The z output structure has the following fields:
%    xbest        point where f(x)=fbest
%    fbest        best function value
%    parent       Vector of index of the parent point
%    len          Vector of length of offset from the parent
%    dim          Vector of dimension of offset from the parent
%    S            Vector with rectangle size indices
%    F            Vector with all function values.
%    s            Row vector of all different rectangle sizes, sorted.
%    f            Row vector of minimum function value for each size 
%    iters        Number of iterations
%    nevals       Number of function evaluations
%    nsizes       Number of sizes
%    ibest        location of minimizing value
%    lengths      list of possible side lengths
%    xmid         middle point of initial rectangle
%    xrange       side lengths of initial rectangle
%    convtype     Possible values:
%                    1: maxevals reached
%                    2: maxcuts reached
%                    3: minimum function value reached
%                    4: relative function value reached
%                    5: minimum size reached
%                    6: maximum iterations reached
%
% Requires MEX file directc

% Implements the DIRECT algorithm from
% S. R. Jones, C. S. Perttunen and B. E. Stuckman,
% "Lipschitzian Optimization Without the Lipschitz Constant,"
% Journal of Optimization Theory and Applications, 79(1993): 157-181.
% Additional convergence options and rectangle size options are added.
%
% Coded by Paul L. Fackler, NCSU
% paul_fackler@ncsu.edu

function [xbest,fbest,convtype,z] = direct(func,bounds,options,varargin)
% Check inputs
if nargin < 2;
   disp('Must pass 2 input arguments');
   return;
end
if isempty(bounds) | size(bounds,2)~=2
   disp('Both lower and upper variable bounds are required');
   return;
end
% Determine option values
if nargin<3, options=[]; end
getopts(options, ...
 'maxflag',   0, ...         % set to 1 for max problems, 0 for min problems
 'sizeconst', 0.5,...        % constant on rectangle size function
 'distopt',  1,...           % 1/0 for distance/volume measure of size
 'oneeach',   0,...          % 1 if only 1 rectangle of each size is split per iteration
 'greedy',    0,...          % if 1, use all sizes greater than the min slope size
 'epsilon',   1e-4,...       % global/local weight parameter 
 'maxit',     50,...         % maximum of iterations
 'maxeval',   10000,...      % maximum # of function evaluations
 'maxcuts',    32,...        % maximum number of side divisions
 'fopt',      -realmax,...   % minimum value of function
 'tol',       0,...          % allowable relative error (if 0, function must reach fopt)
 'minlength', 1e-6,...       % stop if best rectangle has all sides 1ess than this
 'minevals',  1000,...       % but must evaluate at least this many points
 'ShowIters', 0);            % print level
if maxflag, minind=-1; if fopt==-realmax, fopt=realmax; end
else, minind=1; 
end
a      = bounds(:,1);
b      = bounds(:,2);
n      = length(a);     % Problem dimension
xrange = (b - a);       % Lengths of boundaries
xmid   = a + xrange/2;  % center of original search space
fbest  = feval(func, xmid, varargin{:});  % Function value at x
nevals = 1;             % number of function evaluations
% Initialize memory
len      = zeros(1,maxeval);
dim      = zeros(1,maxeval);
parent   = zeros(1,maxeval);
S        = zeros(1,maxeval);
F        = zeros(1,maxeval);
S(1)     = 1;                 % Vector of all rectangle sizes
F(1)     = fbest;      % Vector of all function values 
ibest    = 1;                 % location of the minimal value
% Precompute possible side lengths and rectangle sizes
lengths=1./3.^(0:maxcuts);
maxcuts=min(sum(minlength<lengths)+1,maxcuts);
maxcuts=min(maxcuts,255);
minsize=maxcuts*n+1;
if distopt
  s=(9*n-8*(0:n-1)').^sizeconst*(lengths/6);
  s=[s(:)' n.^sizeconst*lengths(end)/6]/n;
  s(n*maxcuts+2:end)=[];
else
  s=sizeconst.^(0:n*maxcuts);
end
maxs = length(s);
f    = zeros(1,maxs);
f(1) = fbest;
ns   = 1;                  % smallest size value yet evaluated
it   = 0;                  % iteration count  
% Set up to call MEX function
parent=int32(parent);
dim=int16(dim);
S=int16(S);
parameters=[epsilon;fopt;tol;minind];
intarray=[nevals;it;ns;ibest;0; ...
  maxit;ShowIters;n;maxeval;maxs;minevals;minsize;oneeach;greedy];
% Call  MEX function
directc(func,xrange,xmid,parent,len,dim,S,F,s,f,parameters,intarray,varargin{:});
% Extract results information
nevals   = intarray(1);    
it       = intarray(2); 
ns       = intarray(3);        
ibest    = intarray(4);
convtype = intarray(5);
xbest    = xmid;
ii = parent(ibest); i=ibest;
while (ii>0) xbest(dim(i))=xbest(dim(i))+len(i); i=ii; ii=parent(ii); end
fbest = minind*F(ibest);
% SAVE RESULTS
if nargout>2
  z.xbest     = xbest;              % point where f(x)=fbest
  z.fbest     = fbest;              % Best function value
  z.F         = minind*F(1:nevals); % All function values
  z.S         = S(1:nevals);        % All rectangle sizes
  z.parent    = parent(1:nevals);   % All rectangle sizes
  z.len       = len(1:nevals);      % All rectangle sizes
  z.dim       = dim(1:nevals);      % All rectangle sizes
  z.s         = s(1:ns);            % vector of possible rectangle sizes
  z.f         = minind*f(1:ns);     % best function values for each size
  z.iters     = it;                 % Number of iterations
  z.FuncEvals = nevals;             % Number of function evaluations
  z.nsizes    = ns;                 % Number of sizes
  z.ibest     = ibest;              % location of minimizing value
  z.lengths   = lengths;            % list of possible side lengths
  z.xmid      = xmid;               % middle point of initial rectangle
  z.xrange    = xrange;             % side lengths of initial rectangle
  z.convtype  = convtype;           % convergence type (see documentation) 
end
return
