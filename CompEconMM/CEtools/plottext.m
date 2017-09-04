%% PLOTTEXT
%
%  Inserts text in figure using a simplified format
%
%  Usage
%    plottext(x,y,txt,ha,va,fs)
%  Input
%    x     : x coordinate
%    y     : y coordinate
%    txt   : text to be inserted
%    ha    : horizontal allignment (def:left)
%    va    : vertical allignment (def:bottom)
%    fs    : fontsize of bullet (def: 18)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function plottext(x,y,txt,ha,va,fs)

if isempty(x)
  xl = xlim;
  x = xl(1);
end
if isempty(y)
  yl = ylim;
  y = yl(1);
end
if nargin<4||isempty(ha), ha='left'; end
if nargin<5||isempty(va), va='bottom'; end
if nargin<6, fs=18; end
text(x,y,txt,'HorizontalAlignment',ha,'VerticalAlignment',va,'FontSize',fs,'Interpreter','Latex')