% MINCOL Returns a vector containing the minimum element in each row of a matrix
% USAGE
%   [z,ind]=mincol(x);
%
% Same as [z,ind]=min(x,[],2);
% Coded as a MEX function to improve speed performance
% Run:
%   mex mincol.c
% to create the MEX file
function [z,ind]=mincol(x);
if nargout>1
  [z,ind]=min(x,[],2);
else
  z=min(x,[],2);
end