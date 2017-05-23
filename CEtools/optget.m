% OPTGET Utility to get previously set function default values 
% USAGE
%   optvalue=optget(funcname,optname,optvalue);
% INPUTS
%   funcname : name of function
%   optname  : name of option
%   optval   : option value
% OUTPUT
%   optval   : the current value of the option
%
% If the named field is not already defined, it will be set to
% optvalue, but optvalue has no effect if the field has already 
% been set. Use OPTSET to change a previously set field.
%
% optget(funcname) returns the current values of the options structure.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function optvalue = optget(funcname,optname,optvalue)

funcname = lower(funcname);
optvar=[funcname '_options'];
eval(['global ' optvar])       % declare a global variable

if nargin==1                   % return the whole option structure
  optvalue=(eval(optvar));
  return
end  

optname  = lower(optname);
% if structure is empty or the named field does not exist
% set to the value passed
if isempty(eval(optvar)) | ~isfield(eval(optvar),optname)
  eval([optvar '.' optname '=optvalue;']); 
% otherwise return the value in the field
else
  optvalue = eval([optvar '.' optname]);
end

