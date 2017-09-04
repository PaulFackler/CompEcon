%% INDEXGRID
%
% Generates Cartesian product of array indices
%
%  Usage
%    [i1,i2,i3,...] = indexgrid(n1,n2,n3,...)
%  Input
%    n1,n2,... : number of indices in each direction
%  Output
%    i1,i2,... : prod(n).1 grid indices
%  Example
%   [i1,i2] = indexgrid(3,2) returns
%   i1 = 1
%        2
%        3
%        1
%        2
%        3
%   i2 = 1
%        1
%        1
%        2
%        2
%        2

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function varargout = indices3(varargin)

d = length(varargin);
I = cell(1,d);
for j=1:d
  I{j} = (1:varargin{j})';
end
I = gridmake(I);

varargout = cell(1,d);
for j=1:d
  varargout{j} = I(:,j);
end