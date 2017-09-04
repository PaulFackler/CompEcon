%% DEMQUA06 Area Under a Curve
%
% Illustrates Trapezoid and Simpson's rules graphically.

% Preliminary tasks
demosetup(mfilename)

% Function to be integrated
f = @(x) 0.5*x - x.^2 + 2*x.^3;


%% Basic Figure Setup

n = 501;
xmin = -1;
xmax =  1;
x = nodeunif(n,xmin,xmax);
xwid = xmax-xmin;
z = zeros(n,1);

y = f(x);
ymin = min(y);
ymax = max(y);
ywid = ymax-ymin;
ymin = ymin-0.2*ywid;
ymax = ymax+0.1*ywid;


%% Area Under Curve

figure
hold on
for i=1:n
  plot([x(i) x(i)],[ymin+0.02 y(i)],'y-','linewidth',2)
end
plot(x,y,'linewidth',4)
plottext(xmax-0.3,ymax-1.2,'$f$');
plotvdash(xmin,f(xmin))
plotvdash(xmax,f(xmax))
plottext((xmin+xmax)/2,-2,'A','center','middle');
plottext(xmin,ymin-0.1,'$a$','center','top');
plottext(xmax,ymin-0.1,'$b$','center','top');
xlim([xmin-0.2*xwid xmax+0.2*xwid])
ylim([ymin ymax])
set(gca,'xtick',[])
set(gca,'ytick',[])


%% Trapezoid Rule

for nnode=[3 5 9]
  
  xnode = nodeunif(nnode,xmin,xmax);
  ynode = f(xnode);
  for i=1:nnode-1
    xnodetmp = xnode(i:i+1);
    c = [xnodetmp.^0 xnodetmp]\f(xnodetmp);
    j = find(xnodetmp(1)<=x&x<=xnodetmp(2));
    z(j) = [x(j).^0 x(j)]*c;
  end
  
  figure
  hold on
  plot(x,y,'b-','LineWidth',4)
  plot(x,z,'r--','LineWidth',4)
  legend('$f$','$\hat f$','Location','NW')
  for i=1:n
    plot([x(i) x(i)],[ymin+0.02 z(i)],'y-','linewidth',2)
  end
  plot(x,y,'b-','LineWidth',4)
  plot(x,z,'r--','LineWidth',4)
  for i=1:nnode
    plotvdash(xnode(i),f(xnode(i)))
    plotbullet(xnode(i),f(xnode(i)))
  end
  switch nnode
    case 3
      plottext(xnode(1),ymin-0.1,'$x_0=a$','center','top')
      plottext(xnode(2),ymin-0.1,'$x_1  $','center','top')
      plottext(xnode(3),ymin-0.1,'$x_2=b$','center','top')
      title('Trapezoid Rule Two Intervals')
    case 5
      plottext(xnode(1),ymin-0.1,'$x_0=a$','center','top')
      plottext(xnode(2),ymin-0.1,'$x_1  $','center','top')
      plottext(xnode(3),ymin-0.1,'$x_2  $','center','top')
      plottext(xnode(4),ymin-0.1,'$x_3  $','center','top')
      plottext(xnode(5),ymin-0.1,'$x_4=b$','center','top')
      title('Trapezoid Rule Four Intervals')
    case 9
      plottext(xnode(1),ymin-0.1,'$x_0=a$','center','top')
      plottext(xnode(2),ymin-0.1,'$x_1  $','center','top')
      plottext(xnode(3),ymin-0.1,'$x_2  $','center','top')
      plottext(xnode(4),ymin-0.1,'$x_3  $','center','top')
      plottext(xnode(5),ymin-0.1,'$x_4  $','center','top')
      plottext(xnode(6),ymin-0.1,'$x_5  $','center','top')
      plottext(xnode(7),ymin-0.1,'$x_6  $','center','top')
      plottext(xnode(8),ymin-0.1,'$x_7  $','center','top')
      plottext(xnode(9),ymin-0.1,'$x_8=b$','center','top')
      title('Trapezoid Rule Eight Intervals')
  end
  xlim([xmin-0.2*xwid xmax+0.2*xwid])
  ylim([ymin ymax])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  
end


%% Simpson's Rule

for nnode=[3 5 9]
  
  xnode = nodeunif(nnode,xmin,xmax);
  ynode = f(xnode);
  
  for i=1:(nnode-1)/2;
    xnodetmp = xnode(2*(i-1)+1:2*i+1);
    c = [xnodetmp.^0 xnodetmp.^1 xnodetmp.^2]\f(xnodetmp);
    j = find(xnodetmp(1)<=x&x<=xnodetmp(3));
    z(j) = [x(j).^0 x(j).^1 x(j).^2]*c;
  end
  
  figure
  hold on
  plot(x,y,'b-','LineWidth',4)
  plot(x,z,'r--','LineWidth',4)
  legend('$f$','$\hat f$','Location','NW')
  for i=1:n
    plot([x(i) x(i)],[ymin+0.02 z(i)],'y-','linewidth',2)
  end
  plot(x,y,'b-','LineWidth',4)
  plot(x,z,'r--','LineWidth',4)
  for i=1:nnode
    plotvdash(xnode(i),f(xnode(i)))
    plotbullet(xnode(i),f(xnode(i)))
  end
  switch nnode
    case 3
      plottext(xnode(1),ymin-0.1,'$x_0=a$','center','top')
      plottext(xnode(2),ymin-0.1,'$x_1  $','center','top')
      plottext(xnode(3),ymin-0.1,'$x_2=b$','center','top')
      title('Simpson''s Rule Two Intervals')
    case 5
      plottext(xnode(1),ymin-0.1,'$x_0=a$','center','top')
      plottext(xnode(2),ymin-0.1,'$x_1  $','center','top')
      plottext(xnode(3),ymin-0.1,'$x_2  $','center','top')
      plottext(xnode(4),ymin-0.1,'$x_3  $','center','top')
      plottext(xnode(5),ymin-0.1,'$x_4=b$','center','top')
      title('Simpson''s Rule Four Intervals')
    case 9
      plottext(xnode(1),ymin-0.1,'$x_0=a$','center','top')
      plottext(xnode(2),ymin-0.1,'$x_1  $','center','top')
      plottext(xnode(3),ymin-0.1,'$x_2  $','center','top')
      plottext(xnode(4),ymin-0.1,'$x_3  $','center','top')
      plottext(xnode(5),ymin-0.1,'$x_4  $','center','top')
      plottext(xnode(6),ymin-0.1,'$x_5  $','center','top')
      plottext(xnode(7),ymin-0.1,'$x_6  $','center','top')
      plottext(xnode(8),ymin-0.1,'$x_7  $','center','top')
      plottext(xnode(9),ymin-0.1,'$x_8=b$','center','top')
      title('Simpson''s Rule Eight Intervals')
  end
  xlim([xmin-0.2*xwid xmax+0.2*xwid])
  ylim([ymin ymax])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  
end

% SAVE FIGURES
printfigures(mfilename)