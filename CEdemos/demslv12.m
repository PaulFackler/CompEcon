% DEMSLV12 Cournot demonstration
disp('DEMSLV12 Cournot demonstration')
close all

n = 50;
q1 = nodeunif(n,0.1,1.5);
q2 = nodeunif(n,0.1,1.5);

for i=1:n
for j=1:n
   z = cournot([q1(j);q2(i)]);
   z1(i,j) = z(1);
   z2(i,j) = z(2);
end
end

warning off

close all 
figure(1)
axis([0.1 1.5 0.1 1.5])
[c,h] = contour(q1,q2,z1,[0 0]);
set(h,'LineStyle','-','LineWidth',2,'Tag','F1','edgecolor',[0 0 0],'CDataMapping','direct');
hold on
[c,h] = contour(q1,q2,z2,[0 0]);
set(h,'LineStyle','-','LineWidth',2,'Tag','F2','edgecolor',[0 0 0],'CDataMapping','direct');
text(.55,1.4,'\pi_1''>0');
text(.85,1.4,'\pi_1''<0');
text(1.2,.75,'\pi_2''<0');
text(1.2,.55,'\pi_2''>0');
title('Solve Cournot Model via Newton Method')
h=xlabel('q_1');
set(h,'VerticalAlignment','cap')
h=ylabel('q_2');
set(h,'VerticalAlignment','bottom')

q = [0.2;0.2];
optset('newton','maxit',1);
optset('newton','maxsteps',0);
for i=1:10
  qnew = newton('cournot',q);
  plot([q(1) qnew(1)],[q(2) qnew(2)],'-');
  plot(q(1),q(2),'*');
  q = qnew;
end
hold off

figure(2)
axis([0.1 1.5 0.1 1.5])
[c,h] = contour(q1,q2,z1,[0 0]);
set(h,'LineStyle','-','LineWidth',2,'Tag','F1','edgecolor',[0 0 0],'CDataMapping','direct');
hold on
[c,h] = contour(q1,q2,z2,[0 0]);
set(h,'LineStyle','-','LineWidth',2,'Tag','F2','edgecolor',[0 0 0],'CDataMapping','direct');
text(.55,1.4,'\pi_1''>0');
text(.85,1.4,'\pi_1''<0');
text(1.2,.75,'\pi_2''<0');
text(1.2,.55,'\pi_2''>0');
title('Solve Cournot Model via Broyden''s Method')
h=xlabel('q_1');
set(h,'VerticalAlignment','cap')
h=ylabel('q_2');
set(h,'VerticalAlignment','bottom')

q = [0.2;0.2];
optset('broydeni','maxit',1);
optset('broydeni','maxsteps',0);
fjac = fdjac('cournot',q);
fjacinv = inv(fjac);
for i=1:10
  [qnew,fval,fjacinv] = broydeni('cournot',q,fjacinv);
  plot([q(1) qnew(1)],[q(2) qnew(2)],'-');
  plot(q(1),q(2),'*');
  q = qnew;
end
hold off

warning on


prtfigs(mfilename,'Cournot Model Solved Using Newton''s Method',1)
prtfigs(mfilename,'Cournot Model Solved Using Broyden''s Method',2)


optset('newton','defaults')
optset('broydeni','defaults')
