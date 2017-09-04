%% DEMDIF02 Error in Finite Difference Derivatives

function demdif02

% Preliminary tasks
demosetup(mfilename)

n = 21;
c = -nodeunif(n,0,15);
c = sort(c);
h = 10.^c;
x = 1;

% One-sided finite difference derivative
u = x+h;
l = x;
d1 = (exp(u)-exp(l))./(u-l)-exp(x);
d1 = log10(abs(d1));
% d1 = smooth(c,d1,0.2,'rloess');
e1 = log10(eps^(1/2));

% Two-sided finite difference derivative
u = x+h;
l = x-h;
d2 = (exp(u)-exp(l))./(u-l)-exp(x);
d2 = log10(abs(d2));
% d2 = smooth(c,d2,0.2,'rloess');
e2 = log10(eps^(1/3));

% Plot finite difference derivatives
figure
hold on
plot(c,d1,c,d2)
plotvdash(e1,[],'b')
plotvdash(e2,[],'r')
xlim([-15 0])
xlabel('$Log_{10}(h)$')
ylabel('$Log_{10}$ Approximation Error')
text(e1,-12,'$\sqrt{\epsilon}$','fontsize',18)
text(e2,-12,'$\sqrt[3]{\epsilon}$','fontsize',18)
legend('One-Sided','Two-Sided','Location','SW')
title('Error in Numerical Derivatives')

%% SAVE FIGURES
printfigures(mfilename)
