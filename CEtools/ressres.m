% RESSRES Residual function for finding steady state in RE models
% See RESS

% Copyright (c) 1997-2001,  Paul L. Fackler & Mario J. Miranda 
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function r=ressres(theta,func,params,e,n,m,p,explicit);
 if explicit
   se=theta(1:n)';
   ze=theta(n+1:end)';
   xe=feval(func,'x',se,[],ze,[],[],[],params{:});
   g=feval(func,'g',se,xe,ze,e,[],[],params{:});
   h=feval(func,'h',se,xe,ze,e,se,xe,params{:});
   r=[se-g ze-h]';
 elseif 0 % not used
   se=theta(1:n)';
   xe=theta(n+1:n+m)';
   ze=theta(n+m+1:end)';
   g=feval(func,'g',se,xe,ze,e,[],[],params{:});
   f=feval(func,'f',se,xe,ze,[],[],[],params{:});
   h=feval(func,'h',se,xe,ze,e,se,xe,params{:});
   r=[se-g f ze-h]';
 else
   se=theta(1:n)';
   xe=theta(n+1:n+m)';
   ze=feval(func,'h',se,xe,[],e,se,xe,params{:});
   g=feval(func,'g',se,xe,ze,e,[],[],params{:});
   f=feval(func,'f',se,xe,ze,[],[],[],params{:});
   r=[se-g f]';
 end
