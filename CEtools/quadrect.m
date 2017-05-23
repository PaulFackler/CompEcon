% QUADRECT Integrates function on a rectangular region in R^n
% USAGE
%    I = quadrect(f,n,a,b,type,varargin);
% INPUTS
%    f    : function name
%    n    : number of sample points
%    a,b  : d-vectors of left and right endpoints for each dimension
%          (note: the problem dimension, d, is determined by a and b
%    type : type of integration scheme
%           'lege' - Gauss-Legendre
%           'cheb' - Gauss-Chebyshev
%           'trap' - trapezoid rule 
%           'simp' - Simpson rule 
%           'N'    - Neiderreiter equidistributed sequence
%           'W'    - Weyl equidistributed sequence
%           'H'    - Haber  equidistributed sequence
%           'R'    - Monte Carlo
%    varargin: additional function parameters
% OUTPUT
%    I    : integral of function on region [a,b]
%
% USES: qnwlege, qnwcheb, qnwtrap, qnwsimp, qnwequi

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function I = quadrect(f,n,a,b,type,varargin)

switch type
case 'lege'
   [x,w] = qnwlege(n,a,b);
case 'cheb'
   [x,w] = qnwcheb(n,a,b);
case 'trap'
   [x,w] = qnwtrap(n,a,b);
case 'simp'
   [x,w] = qnwsimp(n,a,b);
otherwise
   [x,w] = qnwequi(n,a,b,type);
end
I = w.'*feval(f,x,varargin{:});
