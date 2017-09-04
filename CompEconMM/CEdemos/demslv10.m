%% DEMSLV10 Linear Complementarity Problem
% 
%  Illusttates solutions to 1-dimensional linear complementarity problem.

function demslv10

% Preliminary tasks
demosetup(mfilename)

% Possible Solutions to Complementarity Problem, $f$ Strictly Decreasing
figure
basicsubplot(1,'$f(a)>f(b)>0$', 1.5, 0.5)
plotbullet(1.0,0.5,18,'r')
basicsubplot(2,'$f(a)>0>f(b)$', 0.5,-0.5)
plotbullet(0.5,0.0,18,'r')
basicsubplot(3,'$0>f(a)>f(b)$',-0.5,-1.5)
plotbullet(0.0,-0.5,18,'r')
box off

% Possible Solutions to Complementarity Problem, $f$ Strictly Increasing
figure
basicsubplot(1,'$f(a)<f(b)<0$',-1.5,-0.5)
plotbullet(0.0,-1.5,18,'r')
basicsubplot(2,'$f(a)<0<f(b)$',-0.5, 0.5)
plotbullet(0.0,-0.5,18,'r')
plotbullet(0.5, 0.0,18,'r')
plotbullet(1.0, 0.5,18,'r')
basicsubplot(3,'$0<f(a)<f(b)$', 0.5, 1.5)
plotbullet(1.0,1.5,18,'r')
box off


%% SAVE FIGURES
printfigures(mfilename,1)


function basicsubplot(i,tit,yleft,yright)
fs = 14;  % FontSize
lw = 2;   % LineWidth
subplot(1,3,i)
hold on
axis ([0 1 -2 2],'square','off')
title(tit,'Fontsize',fs)
plot([0;1],[yleft,yright],'LineWidth',5);
plot([0;1],[ 0;0],'k','LineWidth',lw);
plot([0;0],[-2;2],'--k','LineWidth',lw);
plot([1;1],[-2;2],'--k','LineWidth',lw);
plottext(-0.05,0,'0','right','middle',fs)
plottext(0,-2.8,'a','center','bottom',fs)
plottext(1,-2.8,'b','center','bottom',fs)