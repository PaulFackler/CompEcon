% DEMSLV13 Illustrates CP problems
disp('DEMSLV13 Illustrates CP problems')
close all


fs=get(0,'DefaultTextFontSize')*(3/4);

figure(1)
subplot(2,2,1)
hold on
title('f''<0, f(a)<0')
plot([0;1],[-0.5,-1.5],'LineWidth',3);
plot([0;1],[ 0;0],'--','LineWidth',2);
plot([0;0],[-2;2],'LineWidth',2);
plot([1;1],[-2;2],'LineWidth',2);
h=text(0,0,{'0'});
set(h,'VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',fs)
h=text(0,-2,{'a'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1,-2,{'b'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
axis ([0 1 -2 2])
axis square
axis off
hold off

subplot(2,2,2)
hold on
title('f''<0, f(b)>0')
plot([0;1],[1.5,0.5],'LineWidth',3);
plot([0;1],[ 0;0],'--','LineWidth',2);
plot([0;0],[-2;2],'LineWidth',2);
plot([1;1],[-2;2],'LineWidth',2);
h=text(0,0,{'0'});
set(h,'VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',fs)
h=text(0,-2,{'a'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1,-2,{'b'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
axis ([0 1 -2 2])
axis square
axis off
hold off

subplot(2,2,3)
hold on
title('f''<0, f(a)>0>f(b)')
plot([0;1],[0.5,-0.5],'LineWidth',3);
plot([0;1],[ 0;0],'--','LineWidth',2);
plot([0;0],[-2;2],'LineWidth',2);
plot([1;1],[-2;2],'LineWidth',2);h=text(0,0,{'0'});
set(h,'VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',fs)
h=text(0,-2,{'a'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1,-2,{'b'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
axis ([0 1 -2 2])
axis square
axis off
hold off

subplot(2,2,4)
hold on
title('f''>0')
plot([0;1],[-0.5,0.5],'LineWidth',3);
plot([0;1],[ 0;0],'--','LineWidth',2);
plot([0;0],[-2;2],'LineWidth',2);
plot([1;1],[-2;2],'LineWidth',2);
h=text(0,0,{'0'});
set(h,'VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',fs)
h=text(0,-2,{'a'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1,-2,{'b'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
axis ([0 1 -2 2])
axis square
axis off
hold off

prtfigs(mfilename,'Alternative Complementarity Problems',1)
