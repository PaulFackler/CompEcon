% DISCRAND Discrete random variable simulator
% USAGE
%   d = discrand(m,p)
% INPUTS
%   m = number of occurances simulated
%   p = probability vector, where states
%       1,2,...,n occur with probabilities 
%       prob(1),prob(2),...,prob(n)
% OUTPUT
%   d = m by 1 vector of random occurances

% Copyright (c) 1997-2010, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function d = discrand(m,p)
  c = [0;cumsum(p(:))];
  u=rand(m,1);
  try    % use mex function
    d = lookup(c,u,3);
  catch  % use Matlab function
    d=floor(interp1q(c,(1:length(c)+1)',u));
  end
