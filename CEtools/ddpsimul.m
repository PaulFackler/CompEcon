% DDPSIMUL  Monte Carlo simulation of discrete-state/action controlled Markov process
% USAGE
%   spath = ddpsimul(pstar,s,N,x)
% INPUTS
%   pstar   : optimal state transition matrix
%   s       : k by 1 vector of initial states
%   N       : number of simulated time periods
%   x       : optimal controls
% OUTPUT
%   spath   : k by N+1 vector of simulated states

% Copyright (c) 1997-2004, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [spath,xpath] = ddpsimul(pstar,s,N,x)

l = length(size(pstar));
n = size(pstar,2);
u = ones(n,1);
k = length(s);

if l==2 % Infinite Horizon Model
  spath=zeros(k,N+1);
  cp=cumsum(pstar,2); 
  for t=1:N+1 
    spath(:,t) = s;
    if t<=N
      r = rand(k,1); 
      s = 1+sum(r(:,u)>cp(s,:),2);
    end  
  end
else    % Finite Horizon Model
  T = size(pstar,3);
  if N>T,
    warning('Request for simulations beyond the problem time horizon are ignored')
  end
  N = min(N,T);
  spath=zeros(k,N+1);
  for t=1:N+1
    spath(:,t) = s;
    if t<=N
      cp=cumsum(pstar(:,:,t),2);
      r = rand(k,1); 
      s = 1+sum(r(:,u)>cp(s,:),2); 
    end
  end
end

if nargout>1
  xpath=zeros(k,N+1);
  xpath(:)=x(spath(:));
end