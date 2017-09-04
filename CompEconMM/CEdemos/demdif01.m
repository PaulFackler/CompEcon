%% DEMDIF01 - Finite Difference Hessian Evaluation Structure

% Preliminary tasks
demosetup(mfilename)

figure
axis([-2 2 -2 2])
set(gca,'xtick',-2:1:2,'xticklabel','')
set(gca,'ytick',-2:1:2,'yticklabel','')
plottext(-1.0,-2,'$x_1-h_1$','center','top',12)
plottext( 0.0,-2,'$x_1$'    ,'center','top',12)
plottext( 1.0,-2,'$x_1+h_1$','center','top',12)
plottext(-2.3,-1,'$x_2-h_2$','center','middle',12)
plottext(-2.3, 0,'$x_2$'    ,'center','middle',12)
plottext(-2.3, 1,'$x_2+h_2$','center','middle',12)
plottext(-1, 1,'$f^{ - +}$','center','middle',18)
plottext( 0, 1,'$f^{ 0 +}$','center','middle',18)
plottext( 1, 1,'$f^{ + +}$','center','middle',18)
plottext(-1,-1,'$f^{ - -}$','center','middle',18)
plottext( 0,-1,'$f^{ 0 -}$','center','middle',18)
plottext( 1,-1,'$f^{ + -}$','center','middle',18)
plottext(-1, 0,'$f^{ - 0}$','center','middle',18)
plottext( 0, 0,'$f^{ 0 0}$','center','middle',18)
plottext( 1, 0,'$f^{ + 0}$','center','middle',18)
title('Evaluation Points for Finite Difference Hessians','fontsize',16)

%% SAVE FIGURES
printfigures(mfilename)
