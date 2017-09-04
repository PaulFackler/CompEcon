function plotvdash(x,y,color,linewidth)
% Plots vertical dashed line from current x-axis to point (x,y)
yl = ylim;
if isempty(y), y=yl(2); end 
if nargin<3||isempty(color), color='k'; end
if nargin<4||isempty(linewidth), linewidth=2; end
plot([x x],[yl(1) y],':','LineWidth',linewidth,'color',color)