%% DEMSLV07 Linear Complementarity Problem Methods
%
% Solve linear complementarity problem on R^8 using various methods with
% randomly generated data.

% Preliminary tasks
demosetup(mfilename)

% Generate test problem data
n  = 8;
z = randn(n,2)-1;
a = min(z,[],2);
b = max(z,[],2);
q  = randn(n,1);
M  = randn(n,n);
M  = -M'*M;
xinit = randn(n,1);

% Solve by applying Newton method to minmax formulation
tic
optget('lcpsolve','type','minmax');
tic
[xm,zm] = lcpsolve(M,q,a,b,xinit);
tm = toc;

% Solve by applying Newton method to semismooth formulation
optget('lcpsolve','type','smooth');
tic
[xs,zs] = lcpsolve(M,q,a,b,xinit);
ts = toc;

% Solve using Lemke's algorithm
tic
[xl,zl] = lcplemke(M,q,a,b);
tl = toc;

% Print table
fprintf('\n\nHundreds of seconds required to solve randomly generated linear\n')
fprintf('complementarity problem on R^8 using Newton and Lemke methods\n')
fprintf('Algorithm           Time      Norm\n')
fprintf('Newton minmax     %6.2f  %8.0e\n',100*tm,norm(minmax(xm,a,b,zm)))
fprintf('Newton semismooth %6.2f  %8.0e\n',100*ts,norm(minmax(xs,a,b,zs)))
fprintf('Lemke             %6.2f  %8.0e\n',100*tl,norm(minmax(xl,a,b,zl)))