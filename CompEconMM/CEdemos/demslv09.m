%% DEMSLV09 Hard Nonlinear Complementarity Problem with Billup's Function
%
% Solve hard nonlinear complementarity problem on R using semismooth and
% minmax rootfinding methods.  Problem involves Billup's function.  Minmax
% formulation fails; semismooth formulation suceeds.

function demslv09

% Preliminary tasks
demosetup(mfilename)

warning off

% Generate problem test data
xinit = 0;
a = 0;
b = inf;

% Solve by applying Newton method to minmax formulation
optset('ncpsolve','type','minmax')
tic
[xm,zm] = ncpsolve(@billups,a,b,xinit);
tm = toc;

% Solve by applying Newton method to semismooth formulation
optset('ncpsolve','type','ssmooth')
tic
[xs,zs] = ncpsolve(@billups,a,b,xinit);
ts = toc;

% Print table header
fprintf('Hundreds of seconds required to solve hard nonlinear complementarity\n')
fprintf('problem using Newton semismooth and minmax formulations\n')
fprintf('                    Time      Norm         x\n');
fprintf('Newton minmax     %6.2f  %8.0e     %5.2f\n',100*tm,norm(minmax(xm,a,b,zm)),xm)
fprintf('Newton semismooth %6.2f  %8.0e     %5.2f\n',100*ts,norm(minmax(xs,a,b,zs)),xs)

figure

subplot(1,2,1)
hold on
x = nodeunif(500,-0.5,2.5);
plot(x,billupss(x),x,billupsm(x))
legend('Semismooth','Minmax','Location','S')
plothdash([],0)
axis([-0.5 2.5 -1.0 1.5],'square')
title('Difficult NCP')
xlabel('$x$')

subplot(1,2,2)
hold on
x = nodeunif(500,-0.03,0.03);
plot(x,billupss(x),x,billupsm(x))
legend('Semismooth','Minmax','Location','S')
plothdash([],0)
axis([-.03 .03 -.04 .06],'square')
title('Difficult NCP Magnified')
xlabel('$x$')


%% SAVE FIGURES
printfigures(mfilename,1)


%% Billups' function
function [fval,fjac] = billups(x)
fval = 1.01-(1-x).^2;
fjac = 2*(1-x);


%% Minmax reformulation of Billups' function
function [fval,fjac] = billupsm(x)
[fval,fjac] = billups(x);
if nargout<2
  fval = minmax(x,0,inf,fval,fjac);
else
  [fval,fjac] = minmax(x,0,inf,fval,fjac);
end


%% Semismooth reformulation of Billups' function
function [fval,fjac] = billupss(x)
[fval,fjac] = billups(x);
if nargout<2
  fval = ssmooth(x,0,inf,fval,fjac);
else
  [fval,fjac] = ssmooth(x,0,inf,fval,fjac);
end