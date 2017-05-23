% GENETIC finds the global maximum of a function.
%   [x,y,funcevals,results]=genetic(f,x0,options,varargin)
% INPUTS
%   f:          user-supplied function to be maximized
%   x0:         dx2 matrix of bounds [a b]: a<=x<=b
%   varargin:   additional (optional) arguments for f
% OUTPUTS
%   x         : vector of x coordinates in [a,b] that maximizes f(x). 
%   y         : function 'f' evaluated at x.
%   flag      : convergence flag 
%                   0: with tol of opt
%                   1: maxgen exceeded
%                   2: maxgenlast exceeded
%   results   : Field structure with:
%                   x        : all points from last generation
%                   y        : function values from last generation
%                   iters    : number of iterations
%                   gen_last : generation of last improvement
%                   flag     : convergence flag
% OPTIONS  
%   results       results structure from previous run or matrix of initial
%   points
%   maxflag       set to 1 for max problems, 0 for min problems [1]
%   gensize       size of each generation [100]
%   replication   percent of children replicated [10]
%   reproduction  percent of children reproduced [85]
%   degmut        dx1 vector of standard deviations [ones]
%   muttype       type of mutation, 1=persistant, 2=degenerative  [1]
%   inbox         0/1 1 if points must remain inside the bounds [0]
%   exchange      percent of reproduction based on DNA exchange [0]
%   fopt          terminate if within a tolerance of fopt [-inf]
%   tol           termination tolerance [10^-8]
%   numits        number of iterations to perform [200]
%   maxgen        maximum # of generations [200]
%   maxgenlast    maximum # of generations w/o improvement [100]
%   maxtime       maximum time [500]
%   maxtimelast   maximum time w/o improvement [400]
%   showiter      = 1 (0) if you do (not) want to see each generation [1]
%   vecflag       = 1 (0) if f does (not) accept point matrices [0]

% Coded by Nick Kuminoff with modifications by Paul Fackler
% Copyright (c) 2005-2008, Paul L. Fackler & Nick Kuminoff, NCSU
% paul_fackler@ncsu.edu
 
function [x,y,flag,results]=genetic(f,x0,options,varargin)
%--------------------------------------------------------------------------
% (1) SET OPTIONS & PARAMETERS FOR PARENT DISTRIBUTION
%--------------------------------------------------------------------------
if nargin<3, options=[]; end
getopts(options, ...
 'results',      [],...     % results structure or initial set of points
 'maxflag',       0,...     % set to 1 for max problems, 0 for min problems
 'gensize',     100,...     % size of each generation
 'replication',  10,...     % percent of children replicated
 'reproduction', 85,...     % percent of children reproduced
 'exchange',      0,...     % percent of reproduction based on DNA exchange
 'repfunc',      [],...     % specialized reproduction function
 'degmut',        1,...     % degree of mutation
 'muttype',       1,...     % type of mutation
 'mutfunc',      [],...     % specialized mutation function
 'inbox',         0,...     % 1 if points must remain inside the bounds 
 'fopt',       -inf,...     % termination level
 'tol',        1e-8,...     % termination tolerance
 'numits',      200,...     % number of iterations to perform
 'maxgen',      200,...     % maximum # of generations
 'maxgenlast',  100,...     % maximum # of generations w/o improvement
 'maxtime',     500,...     % maximum time 
 'maxtimelast', 400,...     % maximum time w/o improvement
 'showiter',      0,...     % = 1 (0) if you do (not) want to see each iter 
 'vecflag',       0);       % = 1 (0) if f does (not) accept point matrices
if ~isempty(repfunc), existrepfunc=true;
else                  existrepfunc=false;
end
if ~isempty(mutfunc), existmutfunc=true;
else                  existmutfunc=false;
end
a=x0(:,1);
b=x0(:,2);
if any(a>=b)
  error('Bounds on x are inconsistent')
end 
if isempty(results) | isnumeric(results)
  if isempty(results)
    x=diagmult(b-a,rand(size(b,1),gensize))+a(:,ones(1,gensize));
  else
    x=results; clear results;
    gensize=size(x,2);
  end
  if vecflag==0,                                         % loop required  
    for j=1:gensize
      y(j)=feval(f,x(:,j),varargin{:});                  % fill y
    end
  else                                                   % no loop required  
    y=feval(f,x,varargin{:});                            % fill y                
    y=y(:)'; 
  end 
  iters=0;
  funcevals=gensize;
  gen_last_improvement=0; 
  [y,sortind]=sort(y);                               % points ordered by y
  if maxflag==0, 
    sortind=flipud(sortind(:)); 
    y=flipud(y(:)); 
  end                                                % order to minimize
  x=x(:,sortind);                                    % sort x 
else
  x=results.x;
  y=results.y;
  gensize=size(x,2);
  iters=results.iters;
  funcevals=results.funcevals;
  gen_last_improvement=results.gen_last;
  clear results
end
% initial value for f_old
if maxflag==1
  f_old=-inf;                               
else
  f_old=inf;
end   
replication=round(replication*gensize/100); % ensure integer value
rep=floor(reproduction*gensize/100);        % ensure integer value
mut=gensize-replication-rep;                % determine mutation 
if mut<0, error('replication+reproduction must be <=100'), end
parentdraws=2*rep + mut;                    % number of required parents
parents=zeros(length(a),parentdraws);       % create parent matrix
distrib=cumsum(sqrt(1:1:gensize)');         % create ranking
distrib=[0;distrib/distrib(end)];           % CDF of points
% obtains indices of randomly selected parents
getpar=@(u) floor(interp1q(distrib,(1:gensize+1)',u));
if degmut==1,  degmut=ones(length(a),1);   end   % initialize mutation 
if muttype==1, degmut1=std(x,1,2).*degmut; end   % st.dev. mutation
flag=3;
if numits>maxgen-iters, numits=maxgen-iters; end
%--------------------------------------------------------------------------
% MAIN LOOP 
%--------------------------------------------------------------------------
for i=iters+1:iters+numits                     % loop over each generation 
 %-------------------------------------------------------------------------
 % USE PARENT DISTRIBUTION TO PRODUCE NEXT GENERATION        
 %-------------------------------------------------------------------------
 parents=x(:,getpar(rand(parentdraws,1)));
 if existrepfunc   
   x(:,1:rep)=repfunc(parents,rep);               % reproduction
 else
   x(:,1:rep)=reproduce(parents,rep,exchange);               % reproduction
 end
 if existmutfunc   
   x(:,rep+1:rep+mut)=mutfunc(parents,mut); % mutation
 else
   if muttype==2, degmut1=std(x,1,2).*degmut; end            % st.dev. mutation
   x(:,rep+1:rep+mut)=mutate(parents,mut,degmut1,inbox,a,b); % mutation
 end
 %-------------------------------------------------------------------------
 % EVALUATE POINTS AND ORDER THEM
 %-------------------------------------------------------------------------
  if vecflag==0,                                     % loop required  
    for j=1:gensize
      y(j)=feval(f,x(:,j),varargin{:});              % fill y
    end
  else                                               % no loop required  
    y=feval(f,x,varargin{:});                        % fill y                
    y=y(:)'; 
  end      
  funcevals=funcevals+gensize;    
  [y,sortind]=sort(y);                               % points ordered by y
  if maxflag==0, 
    sortind=flipud(sortind(:)); 
    y=flipud(y(:)); 
  end                                                % order to minimize
  x=x(:,sortind);                                    % sort x  
 %-------------------------------------------------------------------------
 % DISPLAY RESULTS
 %-------------------------------------------------------------------------
 f_best=y(gensize);
 if f_best~=f_old, gen_last_improvement=i; end
 f_old=f_best;
 if showiter
    fprintf('Gen: %4i   gen_last_improvement: %4i   f_best: %15.10f   f_avg: %15.10f\n',...
         i,gen_last_improvement,f_best,mean(y)); 
 end
 %-------------------------------------------------------------------------
 % TEST STOPPING CRITERIA  
 %-------------------------------------------------------------------------
 
if i>=maxgen, flag=1; end
 if abs(f_best-fopt)<tol, flag=0; break, end
 if abs(f_best-mean(y))<tol; flag=4; break; end
 if i-gen_last_improvement>=maxgenlast, flag=2; break, end
end
%--------------------------------------------------------------------------
% PREPARE RESULTS
%--------------------------------------------------------------------------
if nargout>3
  results.x=x;                  % All sampled points in last generation
  results.y=y ;                 % All function values
  results.iters=i;
  results.funcevals=funcevals;
  results.gen_last=gen_last_improvement;
  results.flag=flag;            % Reason algorithm stopped
end
x=x(:,gensize);                 % best point
y=y(gensize);                   % f value at best point




% REPRODUCE breeds points to create the next generation.
%   children=reproduce(parents,rep,exchange)
% INPUTS
%   parents:    matrix of parent points
%   rep:        number of children to create  
%   exchange:   percent of reproduction based on DNA exchange      
% OUTPUTS
%   children:   matrix of children points for the next generation                   
% NOTES
%   (1) This function is called by genetic.  
%   (2) This function uses parents in descending order.  reproduce's 
%       sister function, mutate, uses parents in ascending order.
%       Between the two functions, every parent is used once.
function children=reproduce(parents,rep,exchange);   
%--------------------------------------------------------------------------
% (1)  PREPARE PARENTS FOR MATING       
%--------------------------------------------------------------------------
DNAex=round(exchange*rep/100);              % # children with exchanged DNA
hyp=rep-DNAex;                              % # children from hypercube 
%--------------------------------------------------------------------------
% (2)  GENETICALLY ENGINEER DNA EXCHANGE       
%--------------------------------------------------------------------------
temp=reshape(parents(:,1:2*DNAex),[],2);        % prepare parents for mating
swap=round(rand(size(temp,1),1));               % generate random 1-0 vector
kids1=temp(:,1).*swap + temp(:,2).*(1-swap);    % exchange DNA    
kids1=reshape(kids1,[],DNAex);                  % make kids1
%--------------------------------------------------------------------------
% (3)  MAKE HYPERCUBE KIDS       
%--------------------------------------------------------------------------
temp=reshape(parents(:,2*DNAex+1:2*rep),[],2);  % prepare parents for mating 
swap=rand(size(temp,1),1);                      % generate random 1-0 vector
kids2=(sum(temp,2) + diff(temp,1,2).*swap)/2;   % find points on hypercube
kids2=reshape(kids2,[],hyp);                    % make kids2
%--------------------------------------------------------------------------
% (4)  PREPARE OUTPUT       
%--------------------------------------------------------------------------
children=[kids1 kids2];




% MUTATE stimulates mutation in children.
%   children=mutate(parents,mut,degmut,a,b)
% INPUTS
%   parents:    matrix of parent points
%   mut:        number of mutated children to create  
%   degmut:     dx1 vector of standard deviations to create mutation  
%   inbox:      0/1 - 1 is children must be inside [a,b]
%   a,b:        optimal bounds - projects to the           
% OUTPUTS
%   children:   matrix of children points for the next generation                   
% NOTES
%   (1) This function is called by genetic.  
%   (2) This function uses parents in ascending order, starting from the
%       bottom.  mutate's sister function, reproduce, uses 
%       parents in descending order, starting from the top. Between the two 
%       functions, every parent is used once.
function children=mutate(parents,mut,degmut,inbox,a,b);   
%--------------------------------------------------------------------------
% (1)  TAKE DRAWS FROM NORMAL TO CREATE MUTATIONS IN CHILDREN
%--------------------------------------------------------------------------
d=length(degmut);
mutation=diagmult(degmut,randn(d,mut));             % create mutations
children=parents(:,end-mut+1:end)+mutation;         % mutated children
if inbox, children=min(max(children,a*ones(1,mut)),b*ones(1,mut)); end

