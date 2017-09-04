%% DEMMATH04 Standard Copulas
%
% Draws contour plots and scatter diagrams for Clayton and Gumbel copulas.

% Preliminary tasks
demosetup(mfilename)

tau = 0.7;
nc = 5;
zlim = 2.5;
n = 800;
u = linspace(0.01,0.99,n);
[U1,U2] = meshgrid(u,u);
U = [U1(:) U2(:)];

%% Plot Contours

figure

subplot(1,2,1)
hold on
c = copulapdf('Clayton',U,2*tau/(1-tau));
c = reshape(c,n,n);
z = icdf('Normal',u,0,1);
f = pdf('Normal',z,0,1);
c = c.*kron(f',f);
contour(z,z,c,nc,'Linewidth',2)
xlim([-zlim zlim]);
ylim([-zlim zlim]);
set(gca,'XTick',[-2 0 2])
set(gca,'YTick',[-2 0 2])
plothdash([],0)
plotvdash(0,[])
title('Clayton')
xlabel('$z_1$')
ylabel('$z_2$')
axis square

subplot(1,2,2)
hold on
c = copulapdf('Gumbel',U,1/(1-tau));
c = reshape(c,n,n);
z = icdf('Normal',u,0,1);
f = pdf('Normal',z,0,1);
c = c.*kron(f',f);
contour(z,z,c,nc,'Linewidth',2)
xlim([-zlim zlim]);
ylim([-zlim zlim]);
set(gca,'XTick',[-2 0 2])
set(gca,'YTick',[-2 0 2])
plothdash([],0)
plotvdash(0,[])
title('Gumbel')
xlabel('$z_1$')
ylabel('$z_2$')
axis square


%% Scatter Plot

figure
n = 100;

subplot(1,2,1)
hold on
u = copularnd('Clayton',2*tau/(1-tau),n);
z = icdf('Normal',u,0,1);
scatter(z(:,1),z(:,2),'*')
xlim([-zlim zlim]);
ylim([-zlim zlim]);
set(gca,'XTick',[-2 0 2])
set(gca,'YTick',[-2 0 2])
plothdash([],0)
plotvdash(0,[])
xlabel('$\tilde z_1$')
ylabel('$\tilde z_2$')
title('Clayton')
axis square

subplot(1,2,2)
hold on
u = copularnd('Gumbel',1/(1-tau),n);
z = icdf('Normal',u,0,1);
scatter(z(:,1),z(:,2),'*')
xlim([-zlim zlim]);
ylim([-zlim zlim]);
set(gca,'XTick',[-2 0 2])
set(gca,'YTick',[-2 0 2])
plothdash([],0)
plotvdash(0,[])
xlabel('$\tilde z_1$')
ylabel('$\tilde z_2$')
title('Gumbel')
axis square


%% SAVE FIGURES
printfigures(mfilename,1)