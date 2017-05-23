% DEMDIF01 Plot to illustrate finite difference Hessian evaluation 
function demdif01
close all

disp(' ')
disp('DEMDIF01 Plot to illustrate finite difference Hessian evaluation')

figure(1)
h=plot([-1 0 1 -1 0 1 -1 0 1],[-1 -1 -1 0 0 0 1 1 1],'.');
set(h,'markersize',20)
axis([-2 2 -2 2])
set(gca,'xtick',-2:1:2)
set(gca,'xticklabel','')
set(gca,'ytick',-2:1:2)
set(gca,'yticklabel','')

h=text(-1,-2,'x_1-h_1');
set(h,'VerticalAlignment','top','HorizontalAlignment','center')
h=text(0,-2,'x_1');
set(h,'VerticalAlignment','top','HorizontalAlignment','center')
h=text(1,-2,'x_1+h_1');
set(h,'VerticalAlignment','top','HorizontalAlignment','center')

h=text(-2.1,-1,'x_2-h_2');
set(h,'VerticalAlignment','middle','HorizontalAlignment','right')
h=text(-2.1,0,'x_2');
set(h,'VerticalAlignment','middle','HorizontalAlignment','right')
h=text(-2.1,1,'x_2+h_2');
set(h,'VerticalAlignment','middle','HorizontalAlignment','right')

text(-1.075,1.2,'f^{ - +}')
text(-0.075,1.2,'f^{ 0 +}')
text(0.925,1.2,'f^{ + +}')

text(-1.075,-1.2,'f^{ - -}')
text(-0.075,-1.2,'f^{ 0 -}')
text(0.925,-1.2,'f^{ + -}')

text(-1.35,0,'f^{ - 0}')
text(0.075,0,'f^{ 0 0}')
text(1.075,0,'f^{ + 0}')

set(title('Evaluation Points for Finite Difference Hessians'),'fontsize',16)

hh=get(gca,'children');
for i=1:length(hh)
  if strcmp(get(hh(i),'type'),'text')
    set(hh(i),'fontsize',14);
  end
end

prtfigs(mfilename,'Evaluation Points for Finite Difference Hessians',1)