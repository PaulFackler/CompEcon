%% SSMOOTH
%
%  Reformulates an MCP as a semismooth function
%
%  Usage
%    [fnew,Jnew] = ssmooth(x,a,b,f,J)
%  Input
%    x         : n.1 evaluation point
%    a         : n.1 lower bound
%    b         : n.1 upper bound
%    f         : n.1 function value
%    J         : n.n Jacobian (optional, required if Jnew requested)
%  Output
%    fnew      : n.1 transformed function value
%    Jnew      : n.n transformed Jacobian
%  Note
%    The reformulation uses phi-(phi+(f,a-x),b-x)
%  where
%    phi+(y,z)=y+z+sqrt(y.^2+z^2)
%    phi-(y,z)=y+z-sqrt(y.^2+z^2)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [fnew,Jnew] = ssmooth(x,a,b,f,J)

try  % exists('arrayssx','file')
  if nargout<2
    fnew = arrayssx(x,a,b,f);
  else                            % Compute the Jacobian
    [fnew,ff,xx] = arrayssx(x,a,b,f);
    n = size(x,1);
    Jnew = spdiags(ff,0,n,n)*J-spdiags(xx,0,n,n);
  end
catch
  if length(a)==1, a = a+zeros(size(x)); end
  if length(b)==1, b = b+zeros(size(x)); end
  dainf = find(a==-inf);
  dbinf = find(b== inf);
  n = length(x);
  da = a-x;
  db = b-x;
  sq1 = sqrt(f.^2+da.^2);
  pval = f+sq1+da;
  pval(dainf) = f(dainf);
  sq2 = sqrt(pval.^2+db.^2);
  fnew = pval-sq2+db;
  fnew(dbinf) = pval(dbinf);
  if nargout==2
    dpdy = 1+f./sq1;   dpdy(dainf)=1;  % y=f, z=da
    dpdz = 1+da./sq1;   dpdz(dainf)=0;
    dmdy = 1-pval./sq2; dmdy(dbinf)=1;  % y=pval, z=db
    dmdz = 1-db./sq2;   dmdz(dbinf)=0;
    ff = dmdy.*dpdy;                    % ff =  ds/df
    xx = dmdy.*dpdz+dmdz;               % xx = -ds/dx
    Jnew = spdiags(ff,0,n,n)*J-spdiags(xx,0,n,n);
  end
end