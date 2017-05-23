% KERNEL Computes a kernel estimate of a PDF
% USAGE:
%   f=kernel(x,xi,h,ktype);
% INPUTS:
%   x     : an mx1 vector of evaluation points
%   xi    : an nx1 vector of observations
%   h     : bandwidth (optional)
%   ktype : 0 = Gaussian, 1 = Epanechnikov (optional, default=1)
% OUTPUT
%   f     : an mx1 vector of estimates of the PDF
%
% Example:
%   xi=randn(1000,1);
%   x=linspace(-5,5,101)';
%   figure(1); plot(x,[normpdf(x) kernel(x,xi,.4)]);
%
% USES: KERNELX

% Note: KERNELX is a C-mex implementation. If it is not being
% used you should compile it. The C implementation is faster
% and uses less memory than the m-file version.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [f,F]=kernel(x,xi,h,ktype)

if nargin<4 | isempty(ktype)
  ktype = optget('kernel','ktype','epanechnikov');
end

if ~isnumeric(ktype)
  switch lower(ktype(1))
    case 'e', ktype = 1;
    case 'g', ktype = 0;
    otherwise, error('Invalid kernel type')
  end
end

if nargin<3 | isempty(h)
  h=0;
end

if nargout>1
  [f,F] = kernelx(x,xi,h,ktype);
else
  f     = kernelx(x,xi,h,ktype);
end
