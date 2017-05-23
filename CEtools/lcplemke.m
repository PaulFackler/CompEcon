% LCPLEMKE Solves linear complementarity problem using Lemke's algorithm
%     z = M*x + q
%     a <= x <= b
%     x_i > a_i => z_i => 0
%     x_i < b_i => z_i =< 0
%
% USAGE
%   [x,z] = lcplemke(M,q,a,b,x);
% INPUTS
%   M       : n by 1 matrix
%   q       : n by 1 vector
%   a       : n by 1 vector, left bound on x
%   b       : n by 1 vector, right bound on x
%   x       : n by 1 vector, initial guess for solution
% OUTPUTS
%   x       : solution to lcp
%   z       : function value at x
%
% Setable options (use OPTSET):
%   tol     : convergence tolerance
%   maxit   : maximum number of iterations

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,z] = lcplemke(M,q,a,b,x)

maxit    = optget('lcplemke','maxit',100);
tol      = optget('lcplemke','tol',sqrt(eps));

if nargin < 5, x=q; end

MM = [-M eye(size(M)) ; -eye(size(M)) zeros(size(M))];
d = length(x);
Minv = inv(M);
MM = [-Minv Minv ; Minv -Minv];
qq = [-Minv*q-a ; Minv*q+b];
z0 = [ max(-M*x-q,0) ; max(M*x+q,0) ];
z = lemke(MM,qq,z0);
z = z(d+1:2*d)-z(1:d);
x = Minv*(z-q); 
