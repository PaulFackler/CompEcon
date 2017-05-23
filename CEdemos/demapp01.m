% DEMAPP01 Plot basis functions and standard nodes for major approximation schemes
disp('DEMAPP01 Plot basis functions and standard nodes for major approximation schemes')
close all

% Set degree and domain of interpolation
n = 9;
a = 0;
b = 1;

% Construct refined unidorm plotting grid
x = nodeunif(1001,a,b);

% Plot monomial basis functions
figure(1);
for j=1:n
   subplot(3,3,j); plot(x,x.^(j-1),'k','LineWidth',2);
   axis([ 0 1 -0.05 1.05]); set(gca,'Ytick',[0 1])
   set(gca,'Xtick',[0 1]);
   if j>=7 
     set(gca,'XtickLabel',{'0' '1'});
   else
     set(gca,'XtickLabel',[]); 
   end
end
subplot(3,3,2); title('Monomial Basis Functions','FontSize',14);

% Plot Chebychev basis functions and nodes
S     = fundefn('cheb',n,a,b);
xnode = funnode(S);
phi   = funbas(S,x);

figure(2);
for j=1:n
   subplot(3,3,j); plot(x,phi(:,j),'k','LineWidth',2);
   axis([ 0 1 -1.05 1.05]); set(gca,'Ytick',[-1 1])
   set(gca,'Xtick',[0 1]);
   if j>=7 
     set(gca,'XtickLabel',{'0' '1'});
   else
     set(gca,'XtickLabel',[]); 
   end
end
subplot(3,3,2); title('Chebychev Polynomial Basis Functions','FontSize',14);

figure(3)
plot(xnode,zeros(size(xnode)),'*')
title('Chebychev Nodes','FontSize',14);

% Plot linear spline basis functions
S     = fundefn('spli',n,a,b,1);
xnode = funnode(S);
phi   = full(funbas(S,x));

figure(4);
for j=1:n
   subplot(3,3,j); plot(x,phi(:,j),'k','LineWidth',2);
   axis([ 0 1 -0.05 1.05]); set(gca,'Ytick',[0 1])
   set(gca,'Xtick',[0 1]);
   if j>=7 
     set(gca,'XtickLabel',{'0' '1'});
   else
     set(gca,'XtickLabel',[]); 
   end
end
subplot(3,3,2); title('Linear Spline Basis Functions','FontSize',14);

figure(5)
plot(xnode,zeros(size(xnode)),'*')
axis([ 0 1 -0.05 0.05]); set(gca,'Xtick',[0 1]);
set(gca,'YTick',[],'DataAspectRatio',[4 2 2])
title('Linear Spline Standard Nodes','FontSize',14);

% Plot cubic spline basis functions and nodes
S     = fundefn('spli',n,a,b);
xnode = funnode(S);
phi   = full(funbas(S,x));

figure(6);
for j=1:n
   subplot(3,3,j); plot(x,phi(:,j),'k','LineWidth',2);
   axis([ 0 1 -0.05 1.05]); set(gca,'Ytick',[0 1])
   set(gca,'Xtick',[0 1]);
   if j>=7 
     set(gca,'XtickLabel',{'0' '1'});
   else
     set(gca,'XtickLabel',[]); 
   end
end
subplot(3,3,2); title('Cubic Spline Basis Functions','FontSize',14);

figure(7)
plot(xnode,zeros(size(xnode)),'*')
title('Cubic Spline Standard Nodes','FontSize',14);

% Plot finite-difference basis functions and nodes
S     = fundefn('lin',n,a,b);
xnode = funnode(S);
phi   = full(funbas(S,x));

figure(8);
for j=1:n
   subplot(3,3,j); plot(x,phi(:,j),'k','LineWidth',2);
   axis([ 0 1 -0.05 1.05]); set(gca,'Ytick',[0 1])
   set(gca,'Xtick',[0 1]);
   if j>=7 
     set(gca,'XtickLabel',{'0' '1'});
   else
     set(gca,'XtickLabel',[]); 
   end
end
subplot(3,3,2); title('Finite Difference Spline Basis Functions','FontSize',14);

figure(9)
plot(xnode,zeros(size(xnode)),'*')
title('Finite Difference Standard Nodes','FontSize',14);

prtfigs(mfilename,'Chebychev Nodes',3)
prtfigs(mfilename,'Monomial Basis Functions',1)
prtfigs(mfilename,'Chebychev Polynomial Basis Functions',2)
prtfigs(mfilename,'Linear Spline Basis Functions',4)
prtfigs(mfilename,'Cubic Spline Basis Functions',6)