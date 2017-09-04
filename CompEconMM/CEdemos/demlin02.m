%% DEMLIN02 Ill-Conditioning of Vandermonde Matrices

% Preliminary tasks
demosetup(mfilename)

% Compute approximation error and matrix condition number
warning off
n = (6:100)';
errv = zeros(length(n),1);
conv = zeros(length(n),1);
for i=1:length(n)
  v = vander(1:n(i));
  errv(i) = log10(norm(eye(n(i),n(i))-v\v,2));
  conv(i) = log10(cond(v));
end

% Smooth using quadratic function
X = [ones(size(n)) n n.^2];
errv = X*(X\errv);
conv = X*(X\conv);
errv = 10.^errv;
conv = 10.^conv;

% Plot matrix condition numbers and approximation error
figure
semilogy(n,[conv errv])
legend('Condition Numer of $V$','$||I-V^{-1}V||_\infty$')
title('Ill-Conditioning of Vandermonde Matrix $V$')
xlabel('Order of Vandermonde Matrix $V$')
ylabel('Value')

% SAVE FIGURES
printfigures(mfilename)