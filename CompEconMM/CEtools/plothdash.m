function plothdash(x,y,color,linewidth)
% Plots horizontal dashed line from current y-axis to point (x,y)
xl = xlim;
if isempty(x), x=xl(2); end 
if nargin<3||isempty(color), color='k'; end
if nargin<4||isempty(linewidth), linewidth=2; end
plot([xl(1) x],[y y],':','LineWidth',linewidth,'color',color)