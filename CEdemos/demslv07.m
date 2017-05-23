% DEMSLV07 Demonstrates different NCP methods with Billup's function
disp('DEMSLV07 Demonstrates different NCP methods with Billup''s function')
close all

xinit = 0;
a = 0;
b = inf;

warning off;

fprintf('Solution to Billups NCP\n');
fprintf('name                      time    norm   x\n');
for i=1:2
   tic
   switch i
   case 1
      name = 'Newton semismooth     ';
      optset('ncpsolve','type','smooth');
      optset('ncpsolve','maxsteps',50);
      [x,z] = ncpsolve('billups',a,b,xinit);
   case 2
      name = 'Newton minimax        ';
      optset('ncpsolve','type','minmax');
      optset('ncpsolve','maxsteps',50);
      [x,z] = ncpsolve('billups',a,b,xinit);
   end 
   err = norm(minmax(x,a,b,z));
   fprintf('%s%8.2f%10.1e%7.3f\n',name,[toc err  x])
end

warning on;

optset('ncpsolve','defaults');
