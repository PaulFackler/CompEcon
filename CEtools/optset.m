% OPTSET Utility to set function options 
% USAGE
%   optset(funcname,optname,optvalue)
% INPUTS
%   funcname : name of function
%   optname  : name of option
%   optval   : option value
%
% If optname='defaults' the current setting of the options will be
%    cleared. The next time the function is called, the default
%    options will be restored.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function optset(funcname,optname,optvalue)

optvar = [lower(funcname) '_options'];   % name of global variable
optname  = lower(optname);               % name of option field
if strcmp(optname,'defaults')
  eval(['clear global  ' optvar])        % clears global variable
else
  eval(['global  ' optvar])                % declare global variable
  eval([optvar '.' optname '=optvalue;'])  % set specified field
end

