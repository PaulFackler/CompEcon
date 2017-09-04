%% LINDEF
%
%  Defines parameters for linear basis functions and computes standard 
%  breakpoints for linear spline
%
%  Usage
%    [n,a,b,parms] = lindef(breaks,evennum)
%    [n,a,b,parms] = lindef([a,b],num)
%  Input
%    breaks    : breakpoint sequence
%    evennum   : 1 if breakpoints are evenly spaced, 0 otherwise
%    a         : lower bound of approximation interval
%    b         : upper bound of approximation interval
%    num       : number of evenly spaced breakpoints
%  Output
%    n         : number of basis functions (1 plus the polynomial order)
%    a         : lower bound of approximation interval
%    b         : upper bound of approximation interval
%    parms     : cell array containing n, a, b
%  Note: Standard breakpoints are evenly spaced

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [n,a,b,parms] = lindef(breaks,evennum)

if nargin<1 || isempty(breaks)
  error('Breakpoint information must be passed')
end
if nargin<2
  evennum = 0;
end

n=length(breaks);
if n<2
  error('breakpoint sequence must contain at least two elements')
end
if any(diff(breaks)<=0)
  error('breakpoint sequence must be increasing');
end
if evennum~=0
  if length(breaks)==2
    breaks = linspace(breaks(1),breaks(2),evennum)';
  else
    if length(breaks)<2
      error('breakpoint must have at least two elements')
    end
    if any(abs(diff(diff((breaks))))>5e-15*mean(abs(breaks)))
      error('breakpoint sequence is not evenly spaced')
    end
    evennum = length(breaks);
  end
end

n = length(breaks);
a = breaks(1);
b = breaks(end);
parms = {breaks,evennum};