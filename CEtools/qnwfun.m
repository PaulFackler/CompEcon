% QNWFUN Quadrature nodes and weights associated with a function family
% USAGE
%   [x,w]=qnwfun(fspace)
% INPUT
%   fspace : a function definition structure (see FUNDEF)
% OUPUTS
%   x  : quadrature nodes
%   w  : quadrature weights
%
% If f(x) is approximated by funeval(c,fspace,x), then
%    int_a^b f(x)dx 
% is approximated by w'*funeval(c,fspace,x) and
%    int_a^b f(x)g(x)dx 
% is approximated by w'*(funeval(c,fspace,x).*g(x)) and

function [x,w]=qnwfun(fspace)
x=funnode(fspace);
w=full(funbas(fspace,fspace.b,-1)/funbas(fspace,x))';