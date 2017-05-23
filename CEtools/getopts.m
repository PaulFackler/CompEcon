% GETOPTS Returns options values in an options structure
% USAGE
%   [value1,value2,...]=getopts(options,field1,default1,field2,default2,...)
% INPUTS
%   options  : a structure variable
%   field    : a field name
%   default  : a default value 
% OUTPUTS
%   value    : value in the options field (if it exists) or the default value
%
% Variables with the field names will be created in the caller's workspace
% and set to the value in the option variables field (if it exists) or to the 
% default value.
% 
% Example called from a function:
%   getopts(options,'tol',1e-8,'maxits',100);
% where options contains the single field 'tol' with value equal to 1
% The function have two variable defined in the local workspace, tol with a
% value of 1 and maxits with a value of 100.
%
% If options contains a field name not in the list passed to getopts, a
% warning is issued.

function varargout=getopts(options,varargin)
K=fix(nargin/2);
if nargin/2==K
  error('fields and default values must come in pairs')
end
if isa(options,'struct'), optstruct=1; else, optstruct=0; end
varargout=cell(K,1);
k=0;
ii=1;
for i=1:K
  if optstruct & isfield(options,varargin{ii})
    assignin('caller',varargin{ii},getfield(options,varargin{ii}));
    k=k+1;
  else
    assignin('caller',varargin{ii},varargin{ii+1});
  end
  ii=ii+2;
end
  
if optstruct & k~=size(fieldnames(options),1)
  warning('options variable contains improper fields')
end