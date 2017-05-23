% FISHHR residual function for fish harvesting problem
% Used by FISHH

function [e,c]=fishhr(sstar,B,D,f,phil10,nl)
  c=(B+sstar*D)\(sstar*f);
  e=sstar-phil10*c(1:nl);