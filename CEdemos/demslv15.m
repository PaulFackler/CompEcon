% DEMSLV15 Illustrates minmax and semi-smooth reformulations of CP problems
disp('DEMSLV15 Illustrates minmax and semi-smooth reformulations of CP problems')
close all

xmin = -0.3;
xmax =  0.3;
x = nodeunif(500,xmin,xmax);
z = zeros(size(x));

y = -4*x-2*tanh(x);

a = -0.5;
b =  0.5;
ya = a-x;
yb = b-x;

close all
figure(1)

hold on
title('Minimax Reformulation')
ym = minmax(x,a,b,y);
plot(x,ym,'LineWidth',4);
plot(x,z , ':');
plot(x,y ,'LineWidth',2);
plot(x,ya,'--');
plot(x,yb,'--');
text(-0.27,-0.4,{'a-x'})
text( 0.22, 0.4,{'b-x'})
text(-0.18, 1.3,{'f(x)'})
axis ([xmin xmax -1.5 1.5])
set(gca,'Xtick',[ ])
set(gca,'Ytick',[0])
axis square
%axis off
hold off

figure(2)
hold on
title('Semismooth Reformulation')
ys = smooth(x,a,b,y);
plot(x,ys,'LineWidth',4);
plot(x,z , ':');
plot(x,y ,'LineWidth',2);
plot(x,ya,'--');
plot(x,yb,'--');
text(-0.27,-0.4,{'a-x'})
text( 0.22, 0.4,{'b-x'})
text(-0.18, 1.3,{'f(x)'})
axis ([xmin xmax -1.5 1.5])
set(gca,'Xtick',[ ])
set(gca,'Ytick',[0])
axis square
%axis off
hold off

prtfigs(mfilename,'Alternative Reformulations for Complementarity Problems',[1 2])
