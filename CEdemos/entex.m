% ENTEX Solves the firm entry/exit problem
% USAGE
%   [pstar,ci,ca,fspacei,fspacea]=EntExSol(model,ni,na,lo,hi);
% INPUTS
%  model  : a structure variable with model specifications
%  ni, na : number of collocation nodes for the inactive and active regions
%  lo, hi : lower and upper bounds
%            inactive region is approximated on P in P_h*[lo,1]
%              active region is approximated on P in P_l*[1,hi]
% OUTPUTS
%   pstar : [P_l,P_h] optimal switching points
%   ci, ca : coefficients for the value function
%   fspacei, fspacea : familiiy definition structures for the value function
%
% model should contain the following fields:
%   func      : name of the problem definition file (see below)
%   entrycost
%   exitcost
%   productioncost
%   params         : additional parameters to pass to func
% 
% The function definiton file should have the syntax
%   [rho,mu,sigma]=func(P,additional parameters)

function [pstar,ci,ca,fspacei,fspacea]=EntExSol(model,ni,na,lo,hi)

f=model.func;
I=model.entrycost;
E=model.exitcost;
C=model.productioncost;
params=model.params;

% Define nodal values and basis matrices
fspacei=fundef({'cheb',ni-2,lo,1});
yi=funnode(fspacei);
fspacei=fundef({'cheb',ni,lo,1});
fspacea=fundef({'cheb',na-2,1,hi});
ya=funnode(fspacea);
fspacea=fundef({'cheb',na,1,hi});

% Define differential operators and coefficient matrices
phi0i=funbas(fspacei,0);
phi1i=-funbas(fspacei,1);
phi1a=-funbas(fspacea,1);
phiha=funbas(fspacea,hi,2);
Phii=funbasx(fspacei,yi,[0;1;2]);
Phia=funbasx(fspacea,ya,[0;1;2]);

b0=[zeros(ni-1,1);I;E;zeros(na-2,1)-C;0];
b1=[zeros(ni+1,1);ya;0];

R=[ -funbas(fspacei,1,1)       zeros(1,na)    ; 
         zeros(1,ni)      -funbas(fspacea,1,1)];

% Initialize Pstar and call root finding algorithm
pstar=[1;2];   % initial condition
pstar=broyden('entexres',pstar,[],f,params,yi,ya,...
        Phii,Phia,phi0i,phi1i,phi1a,phiha,b0,b1,R,fspacei,fspacea);

% Obtain coefficient vectors
[res,ca,ci]=entexres(pstar,f,params,yi,ya,...
        Phii,Phia,phi0i,phi1i,phi1a,phiha,b0,b1,R,fspacei,fspacea);