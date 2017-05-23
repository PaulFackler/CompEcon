% CHKFIELDS Checks if a variable S is a valid structure with fields F
% USAGE
%   errcode=chkfields(s,f);
% Returns an error code with values:
%   0 : no errors
%   1 : s is not a structure
%   2 : The number of fields in s differs from the number of elements of F
%   3 : The field list of s does not match that of F
%
% NOTE: to check fields interactively use fields(s)

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function errcode=chkfields(s,f);

errcode=0;
if ~isstruct(s) 
  errcode=1;
else
  ff=fieldnames(s); 
  if any(size(ff)~=size(f)) errcode=2;
  elseif ~all(strcmp(ff,f)) errcode=3;
  end
end

