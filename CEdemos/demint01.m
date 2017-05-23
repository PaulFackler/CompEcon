% DEMINT01 Solves inverse demand problem
  disp(' ')
  disp(' DEMINT01 Solves inverse demand problem')
  
p = 0.25;
for i=1:100
  deltap = (.5*p^-.2+.5*p^-.5-2)/(.1*p^-1.2 + .25*p^-1.5);
  p = p + deltap;
  if abs(deltap) < 1.e-8, break, end
end
fprintf(' price: %10.6g\n',p)

