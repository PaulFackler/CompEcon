% FUNNODE Computes default nodes for a family of functions
% USAGE
%   x=funnode(fspace);
% Returns a cell array containing d vectors or an n vector
%   if the problem dimension is 1.
%
% See also: FUNDEF, FUNBAS, FUNFITF, FUNFITXY, FUNEVAL.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function x=funnode(fspace)

% Determine the problem dimension
d=fspace.d;

if d==1                               % a vector returned when d=1
  x=feval([fspace.bastype{1} 'node'],fspace.parms{1}{:});     
else                                  % a cell array returned otherwise
  x=cell(1,d);                        
  for j=1:d
    x{j}=feval([fspace.bastype{j} 'node'], fspace.parms{j}{:});  
  end
end
