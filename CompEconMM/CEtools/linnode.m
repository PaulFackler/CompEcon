%% LINNODE
%
%  Standard nodes for linear spline
%
%  Usage
%    x = linnode(breaks)
%  Input
%    breaks    : n.1 vector of breakpoints
%  Output
%    x         : n.1 vector of nodes
%  See Also: chebase, funnode.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function x = linnode(breaks,evennum)
x = breaks;