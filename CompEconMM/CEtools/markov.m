%% MARKOV 
%
%  Computes invariant distribution of Markov Chain
%
%  Usage
%    [p,f] = markov(Q)
%  Input
%    Q         : n.n Markov transition probability matrix
%                (non-negative matrix with unit row sums)
%  Output
%    p         : n.k matrix of invariant distributions 
%                for each recurrence class
%                p(i,k)=long-run frequency of visits to state i
%                given the system enters recurrence class k
%    f         : n.n matrix of accessibility probabilities
%                f(i,j)= Prob[state j will be reached from state i]
%
% Reference: "Introduction to Stochastic Processes" by Ethan Cinlar
%   Prentice-Hall, 1975. Chapters 5 and 6.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [p,f] = markov(Q)

n = size(Q,1);

% Error Checking to ensure Q is a valid stochastic matrix
if size(Q,2)~=n
  error('Transition matrix is not square');
end
if any(any(Q<0))
  error('Transition matrix contains negative elements');
end
if any(abs(Q*ones(n,1)-1)>1e-12)
  error('Rows of transition matrix do not sum to 1');
end

% Determine accessibility from i to j
f = zeros(n,n);
for j=1:n
  dr = 1;
  r = spones(Q(:,j));            % a vector of ones where Q(i,j)~=0
  while any(dr)
    dr = r;
    r = spones(Q*r+r);
    dr = r-dr;
  end
  f(:,j) = r;
end

% Determine membership in recurrence classes
% Programming note:
%  f(i,:)=1 for states accessible from i
%  f(:,i)'=1 for states from which i is accessible
%  f(:,i)'.*f(i,:)=1 for states communicating with i (two-way accessibility)
%  If the set of communicating states is the same as the set of accessible
%    states, it forms a recurrence class.
ind = zeros(n,n);
numrec = 0;                 % number of recurrence classes
for i=1:n
  if all(ind(i,:)==0)
    j = f(i,:);             % states accessible from i
    if all((f(:,i)'.*j)==j) % are all accessible states communicating states?
      j = find(j);          % members in class with state i
      k = length(j);        % # of members
      if k>0
        numrec=numrec+1;
        ind(j,numrec) = ones(k,1);
      end
    end
  end
end
ind = ind(:,1:numrec);      % ind(i,j)=1 if state i is in class j

% Determine recurrence class invariant probabilities
p = zeros(n,numrec);
for j=1:numrec
  k = find(ind(:,j));             % members in class j
  nk = length(k);                 % # of members
  % solve Qp=p s.t. 1'p=1
  p(k,j) = [ones(1,nk);(speye(nk)-Q(k,k)')]\[1;zeros(nk,1)];
end

% Analyze transients if second output desired
if nargout>1
  if numrec>1 trans = find(sum(ind')==0);
  else trans = find(ind==0);
  end
  numt = length(trans);            % number of transient states
  % Determine transient absorption and reachability probabilities
  if numt>0
    pp = Q(trans,trans);
    b = zeros(numt,n);
    for j=1:numrec
      k = find(ind(:,j));          % members of class j
      nk = length(k);              % # of members
      if nk==1                     % 1-step prob: transient states to class j
        b(:,k) = Q(trans,k);
      else
        b(:,k) = sum(Q(trans,k)')'*ones(1,nk);
      end
    end
    pp = inv(eye(numt)-pp);
    f(trans,:) = pp*b;              % absorption probabilities
    d = diag(pp)';
    pp = pp./d(ones(numt,1),:);
    f(trans,trans) = pp-diag(1./d); % transient reachability probabilities
  end
end