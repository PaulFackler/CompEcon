%% DEMSLV11 Rootfinding Reformulations of Nonlinear Complementarity Problem
%
% Illustrates minmax and semismooth rootfinding reformulations of nonlinear
% complementarity problem.

% Preliminary tasks
demosetup(mfilename)

xmin = -0.3;
xmax =  0.3;
x = nodeunif(500,xmin,xmax);
z = zeros(size(x));

y = -4*x-2*tanh(x);
a = -0.5;
b =  0.5;
ya = a-x;
yb = b-x;

figure
for i=1:2
  subplot(1,2,i)
  hold on
  set(gca,'Xtick',[])
  set(gca,'Ytick',[])
  if i==1
    title('Minimax Reformulation')
    yr = minmax(x,a,b,y);
  else
    title('Semismooth Reformulation')
    yr = ssmooth(x,a,b,y);
  end
  plot(x,yr,'r','LineWidth',5);
  plot(x,y);
  plot(x,z,'k:');
  plot(x,ya,'k--');
  plot(x,yb,'k--');
  plottext(xmin,-0.3,'$a-x$' ,'left'  ,'top',   14)
  plottext(xmax, 0.3,'$b-x$' ,'right' ,'bottom',14)
  plottext(-0.25,1.3,'$f(x)$','center','middle',14)  
  axis ([xmin xmax -1.5 1.5],'square')
  box off
end

%% SAVE FIGURES
printfigures(mfilename,1)