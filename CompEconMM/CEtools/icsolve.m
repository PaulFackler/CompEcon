%% ICSOLVE
%
%  Solves continuous time impulse control model
%
%  Usage
%    [cv,basis,x] = icsolve(model,x,n,type)
%  Input
%    model     : structure variable (defined below)
%    x         : 2.2 initial control values)
%    n         : degree of approximation
%    typ  e    : type of basis functions to use (default='cheb')
%  Output
%    cv        : n.1 solution coefficients
%    basis     : function definition structure
%    x         : 2.2 optimal control values
%  Model Structure Fields
%    func      : model function file name (see below)  
%    params    : additional parameters to pass to model function
%    xindex    : 2.2 defines nature of the control values, see below
%    R         : marginal rewards (2 by 1)
%    F         : fixed costs
%  Notes
%    The problem is defined on [a,b]. 
%   If a is a trigger, the process jumps to A when a is hit
%   If b is a trigger, the process jumps to B when b is hit
%   where a<=A<=B<=b
%   x is a 2 by 2 matrix composed of [a A;b B]
%   The 2 by 2 variable xindex defines the nature of the control x
%    i j xindex(i,j)  Condition                      
%    1 1    0         a is not a trigger
%    1 2    0
%                    a is a trigger and
%    1 1    1            a is a not a choice variable
%    1 1    2            a is a choice variable
%    1 2    0            A is a not a choice variable
%    1 2    1            A is a choice variable
%
%    2 1    0         b is not a trigger
%    2 2    0
%                    b is a trigger and
%    2 1    1            b is a not a choice variable
%    2 1    2            b is a choice variable
%    2 2    0            B is a not a choice variable
%    2 2    1            B is a choice variable
%  Model Function File Format
%    function out=func(flag,s,additional parameters)
%      switch flag
%        case 'f'
%          Return a matrix f representing the reward function
%        case 'mu'
%          Return a matrix mu representing the drift function
%              of the state transition equation
%        case 'sigma'
%          Return a matrix sigma representing the diffusion function
%              of the state transition equation
%        case 'rho'
%          Return a vector representing the (state contingent)
%              discount rates
%        case 'R+'
%          Return the reward associated with trigger s(1) and target s(2) (s(1)<s(2))
%            [R(s(1),s(2)) R_S0(s(1),s(2)) R_S1(s(1),s(2)) R_S0S0(s(1),s(2))+R_S0S1(s(1),s(2))]
%        case 'R-'
%          Return the reward associated with trigger s(1) and target s(2) (s(1)>s(2))
%            [R(s(1),s(2)) R_S0(s(1),s(2)) R_S1(s(1),s(2)) R_S0S0(s(1),s(2))+R_S0S1(s(1),s(2))]
%      end
%   Each of these should be n by 1 where n=# of rows in s

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [cv,basis,x] = icsolve(model,x,n,type)

if ~exist('type','var'), type='cheb'; end

% Unpack model variables
func = model.func;
params = model.params;
xindex = model.xindex;
F = model.F;

% Get basis matrices and nodal points
m = sum(xindex(:,1)>0);
basis = fundefn(type,n-m,0,1);
z = funnode(basis);
basis = fundefn(type,n,0,1);
Phi0  = funbase(basis,z,0);
Phi1  = funbase(basis,z,1);
Phi2  = funbase(basis,z,2);

phi0 = funbase(basis,[0;1],0);
phi1 = funbase(basis,[0;1],1);
phi2 = funbase(basis,[0;1],2);

% Define xindex variable
if F(1)==0; xindex(1,2) = 0; end
if F(2)==0; xindex(2,2) = 0; end
xindex(:,1) = xindex(:,1)-1;

% Call root finding algorithm (broyden) and get value function coefficients
y = x(xindex>0);
y = broyden(@icres,y,x,func,params,xindex,basis,z,Phi0,Phi1,Phi2,phi0,phi1,phi2,F);
x(xindex>0) = y;  x(F==0,2) = x(F==0,1);
[e,cv] = icres(y,x,func,params,xindex,basis,z,Phi0,Phi1,Phi2,phi0,phi1,phi2,F);

% Adjust basis
basis0 = basis;
basis = fundefn(type,n,x(1,1),x(2,1));

    
%% ICRES
%
%  Residual function for impulse control solver
%
%  Usage
%    [e,c] = icres(y,x,func,params,xindex,basis,z,Phi0,Phi1,Phi2,phi0,phi1,phi2,F)
%  Input
%
%  Used by ICSOLVE

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [e,c] = icres(y,x,func,params,xindex,basis,z, ...
  Phi0,Phi1,Phi2,phi0,phi1,phi2,F)

x(xindex>0) = y;
n = size(z,1);
a = x(1,1);
r = x(2,1)-a;
S = r*z+a;

rho   = feval(func,'rho',S,params{:});
mu    = feval(func,'mu',S,params{:});
H     = spdiags(rho,0,n,n)*Phi0-spdiags(mu/r,0,n,n)*Phi1;
sigma = feval(func,'sigma',S,params{:});
if ~isempty(sigma)
  H = H - spdiags(sigma.*sigma/(2*r^2),0,n,n)*Phi2;
end
f = feval(func,'f',S,params{:});

G = []; g = [];
if xindex(1,1)>=0
  if F(1)==0
    R = feval(func,'R+',[x(1) x(1)],params{:});
    H = [H;phi1(1,:)];
    f = [f;r*R(2)];
    if xindex(1,1)>0, G = [G;phi2(1,:)]; g = [g;r.^2*R(4)]; end
  else
    R = feval(func,'R+',x(1,:),params{:});
    H = [H;phi0(1,:)-funbase(basis,(x(1,2)-x(1,1))/r)];
    f = [f;R(1)-F(1)];
    if xindex(1,1)>0, G = [G;phi1(1,:)]; g = [g;r*R(2)]; end
    if xindex(1,2)>0, G = [G;funbase(basis,(x(1,2)-x(1,1))/r,1)]; g = [g;-r*R(3)]; end
  end
end
if xindex(2,1)>=0
  if F(2)==0
    R = feval(func,'R-',[x(2) x(2)],params{:});
    H = [H;phi1(2,:)];
    f = [f;r*R(2)];
    if xindex(2,1)>0, G = [G;phi2(2,:)]; g = [g;r.^2*R(4)]; end
  else
    R = feval(func,'R-',x(2,:),params{:});
    H = [H;phi0(2,:)-funbase(basis,(x(2,2)-x(1,1))/r)];
    f = [f;R(1)-F(2)];
    if xindex(2,1)>0, G = [G;phi1(2,:)]; g = [g;r*R(2)]; end
    if xindex(2,2)>0, G = [G;funbase(basis,(x(2,2)-x(1,1))/r,1)]; g = [g;-r*R(3)]; end
  end
end

c = H\f;
e = G*c-g;