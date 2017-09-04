%% DEMQUA04 Area Under Normal PDF Using Simpson's Rule

% Preliminary tasks
demosetup(mfilename)

% Standard normal probability density function
f = @(x) sqrt(1/(2*pi))*exp(-0.5*x.^2);

% Compute cumulative probability using Simpson's method
z = 1;
[x,w] = qnwsimp(11,0,z);
prob = 0.5 + w'*f(x);

% Nodes for plotting
n = 500; a = -4; b = 4;
x = nodeunif(n,a,b);

% Plot figure
figure
hold on
for i=1:n
  if x(i)<z, plot ([x(i) x(i)],[0 f(x(i))],'y','LineWidth',2), end
end
plot(x,f(x))
plot([a b],[0 0],'k-')
plot([z z],[0 f(z)],'k-','LineWidth',2)
plottext( 1,-0.01,'$z$','center','top',28);
plottext(-2, 0.20,'$\Pr\left(\tilde Z\leq z\right)$','right','middle',22);
annotation('textarrow',[0.33 0.45],[0.52 0.36]);
set(gca,'ytick',[0.0 0.2 0.4])
axis off

% SAVE FIGURES
printfigures(mfilename)