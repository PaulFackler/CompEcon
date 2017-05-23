% DEMSLV08 Compute fixedpoint of using function iteration and Newton method
disp('DEMSLV08 Compute fixedpoint of using function iteration and Newton method')
close all

% Randomly generate starting point
xinit = rand(1);

warning off;

fprintf('\n')
fprintf('      Time required to compute fixpoint \n')
fprintf('      of f(x)=sqrt(x) via different methods starting at x=%4.2f \n',xinit)
fprintf('name                 time       norm conv      x\n');
for i=1:5
   tic
   switch i
   case 1
      name = 'Newton             ';
      optset('newton','maxsteps',0);
      [x,f] = newton('fslv08b',xinit);
      err = norm(f);
   case 2
      name = 'Newton safeguarded ';
      optset('newton','maxsteps',20);
      [x,f] = newton('fslv08b',xinit);
      err = norm(f);
   case 3
      name = 'Broyden            ';
      optget('broyden','maxsteps',0);
      [x,f] = broyden('fslv08b',xinit);
      err = norm(f);
   case 4
      name = 'Broyden safeguarded';
      optset('broyden','maxsteps',20);
      [x,f] = broyden('fslv08b',xinit);
      err = norm(f);
   case 5
      name = 'Function iteration ';
      [x,f] = fixpoint('fslv08a',xinit);
      err = norm(f-x);
  end 
  fprintf('%s%6.2f   %8.1e   %2i %6.2f\n',name,[toc err (err<1.e-7) x])
end


n=25;
xx=zeros(n,3);

x=0.5;
for k=1:n
  x=fslv08a(x);
  xx(k,1)=x;
end

x=0.5;
optset('broyden','maxit',1)
for k=1:n
  [x,fval]=broyden('fslv08b',x);
  xx(k,2)=x;
end

x=0.5;
optset('newton','maxit',1)
for k=1:n
  [x,fval]=newton('fslv08b',x);
  xx(k,3)=x;
end

disp(' ')
disp(' ')
disp('Solution errors by iteration starting at x=1/2.')
disp('Iteration Function         Broyden          Newton')
fprintf('%2i %15.1e %15.1e %15.1e\n',[1:n;abs(1-xx)']);

% uncomment to print in TEX format
%fprintf('%2i & %15.1e & %15.1e & %15.1e \\\\\n',[1:n;1-xx']);

warning on


optset('newton','defaults')
optset('broyden','defaults')
