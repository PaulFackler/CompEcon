% Incomplete script file that illustrates fundefn and funeval in two dimension.


% f   = @(x)  cos(x(:,1))./exp(x(:,2));
% d1  = @(x) -sin(x(:,1))./exp(x(:,2));
% d2  = @(x) -cos(x(:,1))./exp(x(:,2));
% d11 = @(x) -cos(x(:,1))./exp(x(:,2));
% d12 = @(x)  sin(x(:,1))./exp(x(:,2));
% d22 = @(x)  cos(x(:,1))./exp(x(:,2));


% dfit1  = funeval(c,basis,x,[1 0]);
% dfit2  = funeval(c,basis,x,[0 1]);
% error1 = reshape(dfit1-d1(x),nplot);
% figure; surf(xcoord{1},xcoord{2},error1);
% error2 = reshape(dfit2-d2(x),nplot);
% figure; surf(xcoord{1},xcoord{2},error2);


% dfit11 = funeval(c,basis,x,[2 0]);
% dfit22 = funeval(c,basis,x,[0 2]);
% dfit12 = funeval(c,basis,x,[1 1]);
% error11 = reshape(dfit11-d11(x),nplot);
% figure; surf(xcoord{1},xcoord{2},error11);
% error12 = reshape(dfit12-d12(x),nplot);
% figure; surf(xcoord{1},xcoord{2},error12);
% error22 = reshape(dfit22-d22(x),nplot);
% figure; surf(xcoord{1},xcoord{2},error22);