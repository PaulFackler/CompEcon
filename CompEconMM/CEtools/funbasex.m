%% FUNBASEX
%
%  Creates basis structures for function evaluation (Use funbase to obtain 
%  single basis matrices)
%
%  Usage
%    B = funbasex(basis,x,order,bformat)
%  Input
%    basis   : a structure defining a family of basis functions (default: created by fundefn)
%    x       : mxd matrix or 1xd cell array of columns vectors (default: created by funnode)
%    order   : 1.d matrix (default: zeros(1,d)))
%    bformat : a string: 'tensor', 'direct' or 'expanded'
%  Output
%    B : a basis structure (defined below)
%    x : the computed evaluation points if none are passed
%  Defaults for bformat
%    'expanded' if d=1, otherwise
%    'tensor'   if x is a cell array
%    'direct'   if x is a matrix
%  B-Structure
%    vals     : cell array containing basis data (see exception below)
%    format   : 'tensor', 'direct', 'expanded'
%    order    : orders of differentiation ('expanded' format) or smallest 
%               orders of differentiation and number of bases ('tensor' and 'direct' formats)
%  Note
%    Order Determines the # of basis matrices created:
%      for 'tensor' and 'direct' formats order should be:
%           1xd if only a single basis matrix is needed in each dimension
%           2xd specifying the minimum and maximum orders in each dimension
%           kxd listing of all desired basis matrices (only min and max of order used)
%      for 'expanded' format:
%           kxd listing of all desired basis matrices
%    If d=1 the format will be 'expanded' regardless of bformat

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [B,x] = funbasex(basis,x,order,bformat)

if nargin<1 || ~isstruct(basis)
  error('A coefficient structure must be specified')
end
if nargin<2
  x = [];
end
if nargin<3 || isempty(order)
  order = 0;
end
if nargin<4
  bformat = [];
end

% Determine the problem dimension
d = length(basis.n);

% Expand ORDER if it has a single columns
if d>1 && size(order,2)==1
  order = order*ones(1,d);
end

% Initialize basis structure
m = size(order,1);
if m>1
  minorder = min(order)+zeros(1,d);
  numbases = (max(order)-minorder)+1;
else
  minorder = order+zeros(1,d);
  numbases = ones(1,d);
end

B = struct('vals',[],'order',minorder,'format',bformat);
B.vals = cell(max(numbases),d);

if isempty(x)
  x = funnode(basis);
end

if isempty(bformat)
  if isa(x,'cell')
    bformat = 'tensor';
  else
    bformat = 'direct';
  end
end

if d>1
  if ~isa(x,'cell') && strcmp(bformat,'tensor')
    error('Must pass a cell array to form a tensor format basis structure')
  end
  if isa(x,'cell') && strcmp(bformat,'direct')
    % it would be more efficient in this case to
    % use the tensor form to compute the bases and than
    % to use indexing to expand to the direct form
    % x=gridmake(x); % convert to grid for direct form
  end
end

if isa(x,'cell')
  B.format = 'tensor';
else
  B.format = 'direct';
end

% Compute basis matrices
switch B.format
  case 'tensor'
    for j=1:d
      if (m>1)
        orderj = unique(order(:,j));
      else
        orderj = order(1,j);
      end
      if length(orderj)==1
        B.vals{1,j} = feval([basis.bastype{j} 'base'],basis.parms{j}{:},x{j},orderj);
      else
        B.vals(orderj-minorder(j)+1,j) = feval([basis.bastype{j} 'base'],basis.parms{j}{:},x{j},orderj);
      end
    end
  case 'direct'
    for j=1:d
      if (m>1)
        orderj = unique(order(:,j));
      else
        orderj = order(1,j);
      end
      if length(orderj)==1
        B.vals{1,j} = feval([basis.bastype{j} 'base'],basis.parms{j}{:},x(:,j),orderj);
      else
        B.vals(orderj-minorder(j)+1,j) = feval([basis.bastype{j} 'base'],basis.parms{j}{:},x(:,j),orderj);
      end
    end
end

% d=1, switch to expanded format
if size(B.vals,2)==1
  B.format = 'expanded';
  B.order = order;
  B.vals = B.vals(order+(1-min(order)));
  return
end

% Create expanded format
switch bformat
  case 'expanded'
    B = funbconv(B,order,'expanded');
  case 'direct'
    if isa(x,'cell'), B = funbconv(B,order,'direct'); end
end

%% FUNBCONV
%
%  Converts among basis structure formats; conversion can only go "up",
%  i.e., 'tensor' format can be converted to 'direct' or 'expanded' and
%  'direct' can be converted to 'expanded'. The default is to convert to
%  'expanded'. The default order is the function itself (order 0).

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function B = funbconv(b,order,format)

if chkfields(b,{'vals';'order';'format'})
  disp('Basis structure has an invalid set of fields')
  error(' ')
end

d = size(b.order,2);
if nargin<3 || isempty(format), format='expanded'; end
if nargin<2 || isempty(order);
  if strcmp(format,'expanded')
    order = zeros(1,d);
  else
    order = b.order;
  end
end;

[numbas,d1] = size(order);
if d1~=d
  error('ORDER incompatible with basis structure');
end

if any(min(order,[],1)<b.order)% | any(max(order)>b.maxorder)
  error('Order of derivative operators exceeds those of basis');
end

if strcmp(format,'expanded')
  B.vals = cell(numbas,1);
  B.order = order;
  B.format = format;
  
  if strcmp(b.format,'tensor')
    m=1;n=1;
    for j=1:d
      m=m*size(b.vals{1,j},1);
      n=n*size(b.vals{1,j},2);
    end
    for i=1:numbas
      B.vals{i}=b.vals{order(i,d)-b.order(d)+1,d};
      for j=d-1:-1:1
        B.vals{i} = kron(B.vals{i},b.vals{order(i,j)-b.order(j)+1,j});
      end
    end
  elseif strcmp(b.format,'direct')
    n=1;for j=1:d, n = n*size(b.vals{1,j},2); end
    for i=1:numbas
      B.vals{i} = b.vals{order(i,d)-b.order(d)+1,d};
      for j=d-1:-1:1
        B.vals{i} = dprod(B.vals{i},b.vals{order(i,j)-b.order(j)+1,j});
      end
    end
  elseif strcmp(b.format,'expanded')
    %disp('WARNING: Basis is already in expanded form')
    B = b;
  else
    error('Improper basis format')
  end
elseif strcmp(format,'direct') && strcmp(b.format,'tensor')
  B.vals = cell(numbas,d);
  B.order = order;
  B.format = format;
  ind = cell(1,d);
  for j=1:d
    for i=1:size(b.vals,1)
      if ~isempty(b.vals{i,j})
        ind{j} = (1:size(b.vals{i,j},1))';
        break
      end
    end
  end
  [ind{:}] = gridmake(ind);
  for j=1:d
    for i=1:numbas
      if ~isempty(b.vals{i,j})
        B.vals{i,j} = b.vals{i,j}(ind{j},:);
      end
    end
  end
else
  disp('Not implemented for this option')
  error(' ')
end


%% CHKFIELDS
%
%  Checks if a variable S is a valid structure with fields F
%
%  Usage
%    errcode = chkfields(s,f)
%  Output
%    Returns an error code with values
%    0 : no errors
%    1 : s is not a structure
%    2 : The number of fields in s differs from the number of elements of F
%    3 : The field list of s does not match that of F
%
% Note: to check fields interactively use fields(s)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function errcode = chkfields(s,f)

errcode = 0;
if ~isstruct(s)
  errcode = 1;
else
  ff=fieldnames(s);
  if any(size(ff)~=size(f)) errcode = 2;
  elseif ~all(strcmp(ff,f)) errcode = 3;
  end
end