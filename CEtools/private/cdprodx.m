function a=cdprodx(b,c,ind)
% CDPRODX  Iterated direct product of cell array times a matrix 
% USAGE:
%   a=CDPRODX(B,c,ind);
% The direct product of a set of matrices with the same number of rows is
% equivalent to computing row-wise tensor (Kronecker) products:
%   a(i,:) = (B{ind(1)}(i,:) x ... x B{ind(d)}(i,:))*c
% INPUTS:
%   B : a pxq cell array containing matrices of dimension m by n(i)
%   c : a prod(n(ind(i))) by k matrix
%   ind : d-vector of indices for matrices to select from cell array 
%           (default: all matrices in B, 1:p*q)
% OUTPUT:
%   a : an m by k matrix
%
% See Also: DPROD, KRON, CKRON, CKRONX.

% Copyright (c) 1997, 1999 by Paul L. Fackler
% Coded as a C file.

global CompEcon_MEXwarned

if isempty(CompEcon_MEXwarned)
  disp('Warning: You are using the m-file version of a function that is coded as a C Mex file.')
  disp('  Running this M file version may be significantly slower and more memory intensive.')
  disp('  Use MEXALL to create the executable (MEX or DLL) and make sure it is on the MATLAB path.')
  CompEcon_MEXwarned=1;
end

if nargin<3
  if nargin<2
    error('Must pass two parameters');
  else
    ind=1:numel(b);
  end
end

if ~iscell(b)
  a=b*c;
else  
  d=length(ind);
  a=b{ind(d)}; 
  for i=d-1:-1:1 
    a=dprod(b{ind(i)},a);
  end
  a=a*c;
end

