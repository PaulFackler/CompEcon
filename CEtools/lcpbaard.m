% LCPBAARD Uses block Baard method to solve linear complementarity problem
%     z = M*x + q
%     a <= x <= b
%     x_i > a_i => z_i => 0
%     x_i < b_i => z_i =< 0
% USAGE
%   [x,z] = lcpbaard(M,q,a,b,x);
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

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,z] = lcpbaard(M,q,a,b,x)

tol      = optget('lcpbaard','tol',sqrt(eps));
if nargin < 5, x=q; end

n = length(x);
for i=1:n
   coord{i} = [-1;0;1];
end
nbase = 3^n;
visit = zeros(nbase,1);
powr3 = 3.^[n-1:-1:0];
basis = [];
for it=1:nbase
   w = x+M*x+q;
   i = find(w<=a);
   j = find(a<w & w<b);
   k = find(w>=b);
   l = 2*sum(powr3(k))+sum(powr3(j))+1;
   if visit(l)==1
      if isempty(basis), basis = gridmake(coord); end
      l = min(find(visit==0));
      t = basis(l,:)';
      i = find(t==-1);
      j = find(t== 0);
      k = find(t== 1);
   end
   visit(l)=1;
   if isempty(i)&isempty(k)
      x = -M\q;
   else
      x(i) = a(i);
      x(k) = b(k);
      if ~isempty(j)
         x(j) = -M(j,j)\(q(j)+M(j,[i;k])*x([i;k]));
      end
   end
   z = M*x+q;
   if all(z(i)-tol<0)&all(z(k)+tol>0)&all(a<x+tol)&all(x<b+tol), return, end;
end
warning('Failure to converge in lcpbaard');
