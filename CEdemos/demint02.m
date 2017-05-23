% DEMINT02 Solves nonlinear rational expectations agricultural market model
  disp(' ')
  disp(' DEMINT02 Solves nonlinear rational expectations agricultural market model')

[y,w] = qnwnorm(10,1,0.1);
a = 1;
for it=1:100
  aold = a;
  p = 3 - 2*a*y;
  f = w'*max(p,1);
  a = 0.5 + 0.5*f;
  if abs(a-aold)<1.e-8, break, end
end
fprintf('  a: %10.6g\n',a)
fprintf('  f: %10.6g\n',f)
fprintf(' Ep: %10.6g\n',w'*p)

