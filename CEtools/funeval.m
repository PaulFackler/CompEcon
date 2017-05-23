% FUNEVAL Evaluates multivariate functions with linear bases.
% Evaluates a multivariate function 
%   coefficients g and family definition structure gdef
%   at the points defined by B 
%   at order of differentiation defined by order.
% USAGE
%   y=funeval(g,gdef,B,order);
% INPUTS
%   g       : a coefficient matrix
%   gdef    : a family definition structure (described below)
%   [B]     : a basis structure or
%               an mxd matrix or 1xd cell array of vectors 
%                   default: standard nodes computed by FUNNODE
%   [ORDER] : a d-column matrix (a single column is expanded)
%                   default: zeros(1,d))
% 
% gdef is a structured data type containing fields
%   d        : the dimensional of the input space (number of input variables)
%   n        : the number of coefficients for each dimension
%   bastype  : the type of basis used (e.g., CHEB, SPLI, etc)
%   parms    : a 1xd cell array contianing parameters defining the function
%
% B can be a kxd matrix, a 1xd cell array of column vectors
%   of a basis structure (see FUNBAS for description)
% If B is not a basis structure then FUNEVAL calls FUNBAS to form one,
% which is returned as an optional second output argument. In this case
% B contains X values and may have either one of two forms:
%   A 1xd cell array containing d column vectors 
%   A kxd matrix
% The former will return the function evaluated at all grid points defined
% by the vectors; the latter will return the function at the k points 
% corresponding to the rows of B.
%
% ORDER defines the order of the differential operator to be applied.
% Example: 
%   order=[0 1 1] 
% is the cross partial with respect to the 2nd and 3rd arguments. If 
% multiple evaluations are desired, use a row for each.  
% Example:
%   order=[0 0 0;1 0 0;0 1 0;0 0 1]
% returns the function and its three partial derivatives.
% Negative values of ORDER produce integrals (anti-derivatives),
% normalized to equal 0 at the left endpoints of the approximation
% intervals.
%
% FUNEVAL returns a 3-D array with dimensions representing
%    1: K input evaluation points (K depends on X - see below)
%    2: output dimension (equals # of columns in C)
%    3: derivative dimension (specified by ORDER)
%
%   If B is a KxD matrix 
%     or a basis structure of 'direct' or 'expanded' type 
%       with the 'vals' field containing matrices with k rows
%     then the output of FUNEVAL has K rows
%   If B is a D-dimensional cell array of vectors of length Ki 
%     or a basis structure with D basis matrices with row dimensions Ki 
%     then the output of FUNEVAL has prod(Ki) rows
%   
%   The number of elements in the 3rd dimension of the output is typically
%     the number of rows of ORDER.  If ORDER is omitted, only the function
%     itself is evaluated (no derivatives) and the third dimension has 
%     1 element.  The exception to this rule concerns when a basis 
%     structure of the 'expanded' type is passed. The default is then to 
%     evaluate using all of the bases defined in the structure; K is then
%     equal to length(B.vals).
%
% An optional second output is the basis structure used to evaluate the
% the function.
%
% FUNEVAL makes use of 3 subfunctions FUNEVAL1, FUNEVAL2, and FUNEVAL3
% correspoding to the 3 basis types, 'tensor', 'direct' and 'expanded'.
%
% Note: If B is a 'tensor' type structure or a cell array of column 
%  vectors representing grid values, the output of FUNEVAL is rearranged
%  to make it easier to use, especially for plotting. It is returned as
%  a D+1 dimensional array with the first D-dimensions representing the
%  values of the input arguments and the last dimension representing the
%  output arguments. 
% To find the value of the function at the gridpoint 
%   defined by X(i1,i2,...,id) use
%      y(i1,i2,...,id,:) 
% (the last colon is to get all p outputs).
%  To illustrate, suppose
%    x={(0:0.01:1)' (0:0.025:1)'};
%  This is a 1x2 cell array of column vectors.  Suppose g is a 
%  coefficient structure representing a mapping from R^2 to R^2.
%    FUNEVAL(g,x) will return a 3-D array that is 101x26x2.
% For 2-D plots use 
%   mesh(x{:},FUNEVAL(g,x)')
% Note: transpose is needed because mesh wants x{1} in columns,
%   x{2} in rows.
% If the function has multiple outputs use
%   y=FUNEVAL(g,x);
%   mesh(x{:},y(:,:,k))
% If ORDER has multiple rows, an additional dimension is added to 
% the result.
%
% USES: funbasx, funeval1, funeval2, funeval3
%
% See also: FUNDEF, FUNBAS

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [y,B]=funeval(g,gdef,B,order)

if nargin<1 || (~isa(gdef,'struct') && ~isa(B,'struct'))
  error('Either a coefficient or basis structure must be passed'); 
end
if nargin<3; B=[]; end
if nargin<4, order=0; end

if isempty(g)
  error('Coefficient values not initialized')
end

if ~isa(B,'struct')   % B is not a basis structure so construct one
  d=gdef.d;
  if size(B,2)~=d
    error(['x must have d=' num2str(d) ' columns'])
  end
  if size(order,2)==1
      order=order*ones(1,d); 
  end
  if 0  % legacy code - ignore
  if all(order==0) & isa(B,'double') & all(strcmp(gdef.bastype,'cheb')) %#ok<UNRCH>
    if all(isreal(g(:))) & all(isreal(B(:)))
      y=chebeval(g,B,gdef.n,gdef.a,gdef.b);
      return      
    end
  end
  end
  B=funbasx(gdef,B,order);
else
  if size(order,2)==1, order=order*ones(1,size(B.order,2)); end
end
  
if isa(B.vals,'cell') 
  switch B.format                  % determine the evaluator to call
    case 'tensor'  
      y=funeval1(g,B,order);
    case 'direct'
      y=funeval2(g,B,order);
    case 'expanded'
      y=funeval3(g,B,order);
    otherwise 
      error('Basis structure has an invalid ''format'' field');
  end
else
  y=B.vals*g;
end
