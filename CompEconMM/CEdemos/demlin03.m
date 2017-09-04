%% DEMLIN03 Sparse Linear Equations

% Preliminary tasks
demosetup(mfilename)

N = 800;
M = 100;

AA = rand(N,N);
bb = ones(N,1);
for i=1:N
  for j=1:N
    if abs(i-j)>1, AA(i,j) = 0; end
  end
%   AA(i,i) = sum(AA(i,:));
end

n = [20:30:N]';
ratio = ones(length(n),1);

for k=1:length(n)
  A = AA(1:n(k),1:n(k));
  b = bb(1:n(k));
  tic
  for i=1:M
    x = A\b;
  end
  toc1 = toc;
  S = sparse(A);
  tic
  for i=1:M
    x = S\b;
  end
  toc2 = toc;
  ratio(k) = toc2/toc1;
end

X = [ones(size(n)) n n.^2];
ratio = exp(X*(X\log(ratio)));

% Plot effort ratio
figure
plot(n,ratio)
xlabel('$n$'); 
ylabel('Ratio');
title('Ratio of Sparse to Full Solve Time - Tridiagonal Matrix')


%% SAVE FIGURES
printfigures(mfilename)