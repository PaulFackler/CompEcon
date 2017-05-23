% BISECT Uses method of bisection to find roots for 1-D functions
% Method of bisection for finding the root of a function F, which is known
% to lie on the interval [a,b].  If the function has the same sign
% at both endpoints then an error is printed and execution stops.
% This may be commented out and a problem with an unbracketed root will 
% generate a return equal to one of the endpoints.  Either increasing or
% decreasing functions may be used.  The procedure stops when the change 
% X is within TOL of the root.  
% Note: this can be used by functions that return matrix values.
% USAGE
%   x=bisect(f,a,b,P1,P2,...)
% INPUTS
%   f: a function name (either inline or an M-file) 
%      that takes and returns an nx1 vector or mxn matrix
%   a: a matrix of lower bounds
%   b: a matrix of upper bounds
%   P1,P2,...: optional extra parameters to be passed to F
% OUTPUT
%   x: the roots of f, i.e., f(x)=0
%
% Examples:
%   Suppose the function f is defined in a file as
%       function y=f(x,p1,p2)
%       y=p1+p2.*x.^3;
%   To find roots of f on [-1,-2] and [5,6] use
%       x=bisect('f',[-1;-2],[5,6],p1,p2)
%   Alternatively, use
%       f=inline('p1+p2.*x.^3','x','p1','p2');
%       x=bisect(f,[-1;-2],[5,6],p1,p2);
%
% Setable options (use OPTSET):
%   tol         : tolerence criteria; error less than tol*(b-a) (default=1e-4)
%   checks      : 0 if no checks should be run (default=0)
%   mustbracket : 0 if the root need not be bracketed (default=1) 

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function x=bisect(f,a,b,varargin)
% get options
  tol         = optget('bisect','tol',1e-4);
  checks      = optget('bisect','checks',0);
  mustbracket = optget('bisect','mustbracket',1);

% Perform checks
  if checks
    if nargin<3
      error('At least three parameters must be passed to BISECT');
    end
    if size(a)~=size(b)
      error('In BISECT: Lower and upper ranges must be the same size');
    end
    if any(a>b)
      error('Lower bound greater than upper bound');
    end
  end
  sa=sign(feval(f,a,varargin{:}));
  sb=sign(feval(f,b,varargin{:}));
  if any(sa==sb) & mustbracket
    error('In BISECT: root not bracketed') 
  end

% Initializations  
  dx = 0.5*(b - a);
  tol = dx.*tol;
  x = a + dx;
  dx = sb.*dx;

% Iteration loop
  while any(abs(dx)>tol)
    dx = 0.5*dx;
    x = x - sign(feval(f,x,varargin{:})).*dx;
  end

 return