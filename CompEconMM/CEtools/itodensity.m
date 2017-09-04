%% ITODENSITY
%
%  Derives ergodic distribution of one-dimensional Ito process
%
%  Usage
%    [c,Ex] = itodensity(model,basis,cv)
%  Let
%    dx = dimension of control variable x
%    nb = number of basis functions and collocation nodes (prod(n))
%  Input
%    model  : structured array containing model specifications (see below)
%    basis  : one-dimensional basis defined on state space
%    cv     : nb.1 value function approximant basis function coefficients
%             (optional)
%  Output
%    c      : nb.1 ergodic distribution approximant basis function coefficients
%    Ex     : expected value associated with density 
%  Note
%    Forms an approximation on a bounded interval using
%    p(x)=k*exp(int^x(2*mu(z)/sigma^2(z))dz)/sigma^2(x) where k is a
%    constant that makes p(x) integrate to 1. No checks performed to ensure
%    process admits a proper density.
%  Model Structure
%    The structured array "model" contains fields that specify essential
%    components of the model to be solved (default values for optional
%    fields in parentheses):
%      func      : name of function file (required)
%      params    : model parameters required by function file (empty) 
%  Function File
%    If cv is passed, function file should have same format as function
%    file used by SOCSOLVE. Otherwise the model function file should have
%    the following syntax:
%      function out=funcfile(flag,s,additional parameters)
%      switch flag
%      case 'mu'
%        out = 
%      case 'sigma'
%        out = 
%      end

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [c,Ex] = itodensity(model,basis,cv)

% Set model fields to default if nonexistent (see above)
if ~isfield(model,'func'),    error('Missing Function File');   end  
if ~isfield(model,'params'),  error('Missing Parameter List');  end  
if ~isfield(model,'Vss'),     model.Vss = 0;                    end 

% Unpack model structure
func   = model.func;
params = model.params;
Vss    = model.Vss;

% Derive collocation nodes and interpolation matrix
Phi = funbase(basis);
s = funnode(basis);
ns = size(s,1);

if nargin<3                   % uncontrolled problem
  mu    = feval(func,'mu',s,params{:});
  sigma = feval(func,'sigma',s,params{:});
else                          % controlled problem
  ds = size(s,2);
  Vs = funeval(cv,basis,s,eye(ds));
  if Vss
    Vss = zeros(ns,ds,ds);
    for is=1:ds
      for js=1:ds
        order = zeros(1,ds);
        order(is) = order(is)+1;
        order(js) = order(js)+1;
        Vss(:,is,js) = funeval(cv,basis,s,order);
      end
    end
  else
    Vss = [];
  end
  x     = feval(func,'x',s,[],Vs,Vss,params{:});
  mu    = feval(func,'mu',s,x,[],[],params{:});
  sigma = feval(func,'sigma',s,x,[],[],params{:});
end
sigma = sigma.*sigma;

% Fit ergodic distribution approximation to mu/sigma^2 and integrate
c = Phi\(mu./sigma);
temp = 2*funeval(c,basis,s,-1);
temp = temp-max(temp);                % normalize to avoid overflow

% Fit egrodic distribution approximation to kernel (p)
p = exp(temp)./sigma;
c = Phi\p;
temp = funeval(c,basis,basis.b,-1); % determine constant
c = c./temp;

% Compute expected value (if requested)
if nargout>1
  p = p/temp;
  cc = Phi\(s.*p);
  Ex = funeval(cc,basis,basis.b,-1);
end