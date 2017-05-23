% DEMOPT05 Step length determination
function demopt05
close all

disp(' ')
disp('DEMOPT05 Step length determination')

x=[0.975;.95];

[fx,g0]=banana(x);
d=g0;
optset('optstep','bhhhcone',0.25);
[s1,fx1,errcode,iters1]=optstep(2,'banana',x,fx,g0,d,25);
[s2,fx2,errcode,iters2]=optstep(3,'banana',x,fx,g0,d,25);
[s3,fx3,errcode,iters3]=optstep(4,'banana',x,fx,g0,d,25);
optset('optstep','defaults');
dg=d'*g0;
delta=0.25;
%s=(-0.0010:0.0001:.0025);
s=(-0.0005:0.0001:.0020);
n=length(s);
fs=banana(x(:,ones(n,1))+g0*s);

close all
figure(1)
plot(s,fs-fx,'linewidth',1.5);
hold on
plot(s1,banana(x+s1*d)-fx,'*')
plot(s2,banana(x+s2*d)-fx,'*')
plot(s3,banana(x+s3*d)-fx,'*')
plot(s,s*dg)
ylims=get(gca,'ylim');
plot(s,0,':')
plot([0;0],ylims,':')
s=s(find(s>=0));
plot(s,s*dg*delta,'--')
plot(s,s*dg*(1-delta),'--')
h=text(1.5e-3,ylims(2)/2+4e-5, ...
  {['BHHHSTEP:  s = ' sprintf('%10.8f',s1)];...
   ['  STEPBT:  s = ' sprintf('%10.8f',s2)];...
   ['GOLDSTEP:  s = ' sprintf('%10.8f',s3)]});
fs=get(0,'DefaultTextFontSize')*(3/4);
set(h,'HorizontalAlignment','right','FontSize',fs)
hold off
title('Step Length Determination')
xlabel('s')
ylabel('f(x+sd)')

s={'Iterations';
   ['  BTSTEP: ' num2str(iters2)];...
   ['STEPBHHH: ' num2str(iters1)];...
   ['STEPGOLD: ' num2str(iters3)]};
disp(s)

prtfigs(mfilename,'Step Length Determination',1)
