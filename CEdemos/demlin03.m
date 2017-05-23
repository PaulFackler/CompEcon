% DEMLIN03 Compares various methods for solving sparse linear equations
disp('DEMLIN03 Compares various methods for solving sparse linear equations')
close all

n = 2000;

% Randomly generate problem data
% Matrix A generated in sparse format
c = randn(n,1);
d = randn(n,1);
e = randn(n,1);
d = abs(c)+abs(d)+abs(e);
A = gallery('tridiag',c(2:n),d,e(1:n-1));
b = randn(n,1);

% Matlab Operator (Sparse)
tic
x = A\b;
t1 = toc;

% Gauss Jacobi Iteration (Sparse)
tic
x = gjacobi(A,b);
t2 = toc;

% Gauss Seidel Iteration (Sparse)
tic
x = gseidel(A,b);
t3 = toc;

% Convert A to full format
A = full(A);

% Matlab Operator (Full)
tic
x = A\b;
t4 = toc;

% Gauss Jacobi Iteration (Full)
tic
x = gjacobi(A,b);
t5 = toc;

% Gauss Seidel Iteration (Full)
tic
x = gseidel(A,b);
t6 = toc;

% Print output
fprintf('\n')
fprintf('time needed to solve Ax=b where\n')
fprintf('A is random, diagonally dominant,\n')
fprintf('tridiagonal 100 by 100 matrix.\n')
fprintf('\n')
fprintf('Method                Time\n\n')
fprintf('A\\b Sparse       %9.2f\n',t1)
fprintf('Jacobi Sparse    %9.2f\n',t2)
fprintf('Seidel Sparse    %9.2f\n',t3)
fprintf('A\\b Full         %9.2f\n',t4)
fprintf('Jacobi Full      %9.2f\n',t5)
fprintf('Seidel Full      %9.2f\n',t6)
