%% DPROD    
%
%  Computes direct product of two matrices.
%
%  Usage
%    c = dprod(a,b)
%  Note
%    The direct sum of two matrices with the same number of rows is same as
%    the row-wise tensor (Kronecker) products.

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu%

function c = dprod(a,b)

n = size(a,2); 
p = size(b,2);
c = a(:,ones(p,1)*(1:n)).*b(:,(1:p)'*ones(1,n));