% DEMSLV03 Computes fixedpoint y=[x1^2+x2^3;x1*x2-0.5] via function iteration and Newton methods
disp('DEMSLV03 Computes fixedpoint y=[x1^2+x2^3;x1*x2-0.5] via function iteration and Newton methods')
close all

% Randomly generate starting point
xinit = randn(2,1);

warning off;

fprintf('\n')
fprintf('      Number of flops and time required to compute fixpoint \n')
fprintf('      of f(x1,x2)=[x1^2+x2^3;x1*x2-0.5] via different methods\n')
fprintf('      starting at x1=%4.2f x2=%4.2f\n\n',xinit)
fprintf('name                    time       norm\n');
for i=1:5
   tic
   switch i
   case 1
      name = 'Newton             ';
      optset('newton','maxsteps',0);
      [x,f] = newton('fslv03b',xinit);
      err = norm(f);
   case 2
      name = 'Newton safeguarded ';
      optset('newton','maxsteps',30);
      [x,f] = newton('fslv03b',xinit);
      err = norm(f);
   case 3
      name = 'Broyden            ';
      optset('broyden','maxsteps',0);
      [x,f] = broyden('fslv03b',xinit);
      err = norm(f);
   case 4
      name = 'Broyden safeguarded';
      optset('broyden','maxsteps',30);
      [x,f] = broyden('fslv03b',xinit);
      err = norm(f);
   case 5
      name = 'Function iteration ';
      [x,f] = fixpoint('fslv03a',xinit);
      err = norm(f-x);
  end 
  fprintf('%s%10.4f    %8.1e   %6.2f   %6.2f\n',name,[toc err  x(1) x(2)])
end

warning on;

optset('newton','defaults')
optset('broyden','defaults')