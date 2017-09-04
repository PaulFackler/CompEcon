function plotddash(color,linewidth)
% Plots diagonal dashed line on current axis
xl = xlim;
yl = ylim;
if nargin<3||isempty(color), color='k'; end
if nargin<4||isempty(linewidth), linewidth=2; end
plot([xl(1) xl(2)],[yl(1) yl(2)],':','LineWidth',linewidth,'color',color)