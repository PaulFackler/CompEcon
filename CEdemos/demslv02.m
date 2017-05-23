% DEMSLV02 Computes root of function using Newton, quasi-Newton and Bisection
disp('DEMSLV02 Computes root of function using Newton, quasi-Newton and Bisection')
close all

% Randomly generate starting point
b =  abs(10*randn);
a = -abs(10*randn);
xinit = (a+b)/2;

% Newton's method
tic
for i=1:100
x = newton('fslv02',xinit);
end
time1 = toc;

% Broyden's method
tic
for i=1:100
x = broyden('fslv02',xinit);
end
time2 = toc;

% Bisection method
tic
for i=1:100
x = bisect('fslv02',a,b);
end
time3 = toc;

% Print output
fprintf('\n')
fprintf('Time required to compute roots of exp(-x)-1,\n')
fprintf('100 times using alternative rootfinding methods, \n')
fprintf('starting at a=%5.2f b=%5.2f x=%5.2f\n\n',a,b,xinit)
fprintf('Method             time\n\n')
fprintf('Newton          %5.2f\n',time1)
fprintf('Broyden         %5.2f\n',time2)
fprintf('Bisection       %5.2f\n',time3)