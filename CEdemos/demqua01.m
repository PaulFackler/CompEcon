% DEMQUA01 Plots equi-distributed Sequences in 2-D
disp('DEMQUA01 Plots equi-distributed sequences in 2-D')
close all

type={'Neiderreiter';'Weyl';'Haber';'Random'};
n=2^14;
n=4000;
if 1==1
for i=1:4
  figure(i)
  x=qnwequi(n,[0 0],[1 1],type{i});
  plot(x(:,1),x(:,2),'.','MarkerSize',3);
  title(['2-D ' type{i} ' Type Sequence'])
  xlabel('x_1')
  ylabel('x_2')
  axis square
  set(gca,'xtick',[0 1],'ytick',[0 1])
end
end

if 0
figure(1)
x=qnwequi(n,[0 0],[1 1],'w');
h=plot(x(1,1),x(1,2),'.','erasemode','none');
axis([0 1 0 1])
for i=2:n
  set(h,'Xdata',x(1:i,1),'Ydata',x(1:i,2))
  drawnow
  %if rem(i,100)==0, pause(0.001); end
end
end

n=1000;
for i=1:4
  figure(4+i)
  x=qnwequi(n,[0 0],[1 1],'N');
  plot(x(:,1),x(:,2),'.','MarkerSize',3);
  title(['Neiderreiter Sequence with n=' num2str(n)])
  xlabel('x_1')
  ylabel('x_2')
  axis square
  set(gca,'xtick',[0 1],'ytick',[0 1])
  n=n*2;
end

clear x

prtfigs(mfilename,'Alternative Equidistributed Sequences',[1 2 3 4])
