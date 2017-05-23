% ITOSIMUL  Monte Carlo simulation of a (possibly) controlled Ito process
% USAGE
% For uncontrolled problems:
%   st = itosimul(model,s0,T,n);
% For controlled problems:
%   st = itosimul(model,s0,T,n,cv,fspace);
% INPUTS
%   model   : name of model structure
%   s0      : k by d vector of initial states
%   T       : time horizon (t=linspace(0,T,n))
%   n       : number of time steps
%   cv      : value function coefficients (obtain from SCSOLVE)
%   fspace  : structure for approximating family for value function
% OUTPUT
%   st      : k by d by n matrix of simulated states
% where:
%    k = # of replications
%    d = dimension of the state variable
%    n = # of time steps 
%
% The model structure should contain the following fields:
%   func  : the name of a function file defining the model
%   param : model parameters to be passed to func
%   a,b   : lower and upper bounds on the process
% For uncontrolled problems, the function file should have the syntax
%    out1=funcfile(flag,s,optional additional parameters)
% The values of flag should be
%    'mu'    : the process drift term
%    'sigma' : the process diffusion term
% For controlled problems, the function file should have the syntax
%    out=funcfile(flag,s,x,Vs,optional additional parameters)
% The values of flag should be
%    'x'     : for controlled processes, the optimal control
%                 (this solves f_x+mu_xV_s=0)
%    'g'     : the process drift term
%    'sigma' : the process diffusion term

% Note: this uses an Euler approximation scheme which may not be accurate
%       for large time steps. If any values exceed the bounds([a,b]) those
%       values are reset to the boundary.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function s = itosimul(model,s0,T,n,cv,fspace)

 func=model.func;
 params=model.params;
 a=model.a;
 b=model.b;
 d=size(s0,2); 

 dt=T/n;
 sqdt=sqrt(dt);
 [k,d] = size(s0);
 s = zeros(k,d,n);
 s(:,:,1)=s0;
 st=s0;
 t=0;
 a=repmat(reshape(a,1,d),k,1);      % expand to a kxd matrix
 b=repmat(reshape(b,1,d),k,1);      % expand to a kxd matrix

 if nargin<5             % uncontrolled problem
   for i=2:n
     g = feval(func,'mu',st,params{:});
     sigma = feval(func,'sigma',st,params{:});
     if d==1
       e=sqdt*sigma.*randn(k,1);      % faster for 1-D problems
     else
       e=randn(k,d);
       e=arraymult(sqdt*sigma,e,k,d,d,1);
     end
     st=min(max(real(st+g*dt+e),a),b);
     s(:,:,i) = st;
   end
 else                    % controlled problem
   for i=2:n
     vs=funeval(cv,fspace,st,eye(d));
     x = feval(func,'x',st,[],vs,params{:});
     mu = feval(func,'g',st,x,[],params{:});
     st=st+mu*dt;
     sigma = feval(func,'sigma',st,x,[],params{:});
     if ~isempty(sigma)
       e=randn(k,d);
       if d==1
         e=sqdt*sigma.*e;      % faster for 1-D problems
       else
         e=arraymult(sqdt*sigma,e,k,d,d,1);
       end
       st=st+real(e);
     end
     st=min(max(st,a),b);    s(:,:,i) = st;
   end
 end
