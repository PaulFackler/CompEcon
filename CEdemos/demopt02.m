% DEMOPT02 Displays changes in a simplex
%    creates plot to describe Nelder-Mead
function demopt02
close all

disp(' ')
disp('DEMOPT02 Displays changes in a simplex')

S=[[0.5;0.5] [.65;1.75] [1.5;0.35]];
xy=[0.5 2.25 0.25 2.25];
xy=[0.5 2.225 0.35 2.15];


color0=[.7 .7 .7];
color1=[.2 .2 .2];
rvector=ones(3,1)*(2/2);
rfactor=1+2/2;

x1=S*rvector-rfactor*S(:,1);        % reflection
S1=S; S1(:,1)=x1;
x2=1.5*x1-0.5*S(:,1);                     % expansion
S2=S; S2(:,1)=x2;
x3=0.75*S(:,1)+0.25*x1;                     % contraction
S3=S; S3(:,1)=x3;
S4=S/2+S(:,2)*(ones(1,3)/2);        % shrinkage

fs=get(0,'DefaultTextFontSize')*(3/4);

close all
figure(1)
subplot(2,2,1)
patch(S(1,:),S(2,:),color0)
hold on
patch(S1(1,:),S1(2,:),color1)
hold off
axis(xy)
axis off
h=text(.35,.4,'A');
set(h,'FontSize',fs)
h=text(.55,1.9,'B');
set(h,'FontSize',fs)
h=text(1.5,.2,'C');
set(h,'FontSize',fs)
title('Reflection');

subplot(2,2,2)
patch(S(1,:),S(2,:),color0)
hold on
patch(S2(1,:),S2(2,:),color1)
hold off
axis(xy)
axis off
h=text(.35,.4,'A');
set(h,'FontSize',fs)
h=text(.55,1.9,'B');
set(h,'FontSize',fs)
h=text(1.5,.2,'C');
set(h,'FontSize',fs)
title('Expansion');

subplot(2,2,3)
patch(S(1,:),S(2,:),color0)
hold on
patch(S3(1,:),S3(2,:),color1)
hold off
axis(xy)
axis off
h=text(.35,.4,'A');
set(h,'FontSize',fs)
h=text(.55,1.9,'B');
set(h,'FontSize',fs)
h=text(1.5,.2,'C');
set(h,'FontSize',fs)
title('Contraction');

subplot(2,2,4)
patch(S(1,:),S(2,:),color0)
hold on
patch(S4(1,:),S4(2,:),color1)
hold off
axis(xy)
axis off
h=text(.35,.4,'A');
set(h,'FontSize',fs)
h=text(.55,1.9,'B');
set(h,'FontSize',fs)
h=text(1.5,.2,'C');
set(h,'FontSize',fs)
title('Shrinkage');


%title('Simplex Transformations in the Nelder-Mead Algorithm')

prtfigs(mfilename, 'Simplex Transformations in the Nelder-Mead Algorithm',1)
