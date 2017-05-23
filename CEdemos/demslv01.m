% DEMSLV01 Compute root of Rosencrantz function via Newton and Broyden methods
disp('DEMSLV01 Compute root of Rosencrantz function via Newton and Broyden methods')
close all

optset('newton','defaults')
optset('broyden','defaults')
optset('broydenx','defaults')

% Randomly generate starting point
xinit = [0.71;0.61];

% Newton's method, no backstepping
optset('newton','maxsteps',0);
tic;
x1 = newton('fslv01',xinit);
time1 = toc;

% Newton's method with backstepping
optset('newton','maxsteps',20);
tic;
x2 = newton('fslv01',xinit);
time2 = toc;

% Broyden's inverse method with backstepping
optset('broyden','maxsteps',20);
tic;
x3 = broyden('fslv01',xinit);
time3 = toc;

% Broyden's direct method with backstepping
optset('broydenx','maxsteps',20);
tic;
x4 = broydenx('fslv01',xinit);
time4 = toc;

% Print output
fprintf('\n')
fprintf('Time required to compute roots of Rosencrantz  \n')
fprintf('function via Newton and Broyden methods, \n')
fprintf('starting at x1=%4.2f x2=%4.2f\n\n',xinit)
fprintf('Method                    time   Norm of f     x1     x2\n\n')
fprintf('Newton, no backstep   %8.2f   %8.2e  %5.2f  %5.2f\n',time1,norm(fslv01(x1)),x1)
fprintf('Newton                %8.2f   %8.2e  %5.2f  %5.2f\n',time2,norm(fslv01(x2)),x2)
fprintf('Broyden inverse       %8.2f   %8.2e  %5.2f  %5.2f\n',time3,norm(fslv01(x3)),x3)
fprintf('Broyden direct        %8.2f   %8.2e  %5.2f  %5.2f\n',time4,norm(fslv01(x4)),x4)


optset('newton','defaults')
optset('broyden','defaults')
optset('broydenx','defaults')
