function [x,flag,relres,iter,resvec] = sqmr(A,b,tol,maxit,M1,M2,x0,varargin)
%SQMR    Simplified (Symmetric) Quasi-Minimal Residual Method
%   X = SQMR(A,B) attempts to solve the system of linear equations A*X=B
%   for X.  The N-by-N coefficient matrix A must be square and symmetric.
%   The right hand
%   side column vector B must have length N.  A may be a function AFUN such
%   that AFUN(X) returns A*X and AFUN(X,'transp') returns A'*X.
%
%   SQMR(A,B,TOL) specifies the tolerance of the method.  If TOL is []
%   then SQMR uses the default, 1e-6.
%
%   SQMR(A,B,TOL,MAXIT) specifies the maximum number of iterations.  If
%   MAXIT is [] then QMR uses the default, min(N,20).
%
%   SQMR(A,B,TOL,MAXIT,M) and SQMR(A,B,TOL,MAXIT,M1,M2) use preconditioners M or
%   M=M1*M2 and effectively solve the system inv(M)*A*X = inv(M)*B for X.
%   If M is [] then a preconditioner is not applied.  M may be a function MFUN
%   such that MFUN(X) returns M\X and MFUN(X,'transp') returns M'\X.
%
%   SQMR(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess.  If X0 is []
%   then SQMR uses the default, an all zero vector.
%
%   SQMR(AFUN,B,TOL,MAXIT,M1FUN,M2FUN,X0,P1,P2,...) passes parameters P1,P2,...
%   to AFUN: AFUN(X,P1,P2,...) and AFUN(X,P1,P2,...,'transp') and
%   similarly to the preconditioner functions M1FUN and M2FUN.
%
%   [X,FLAG] = SQMR(A,B,TOL,MAXIT,M1,M2,X0) also returns a convergence FLAG:
%    0 SQMR converged to the desired tolerance TOL within MAXIT iterations.
%    1 SQMR iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 SQMR stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during QMR became too
%      small or too large to continue computing.
%
%   [X,FLAG,RELRES] = SQMR(A,B,TOL,MAXIT,M1,M2,X0) also returns the
%   relative residual NORM(B-A*X)/NORM(B).  If FLAG is 0, RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = SQMR(A,B,TOL,MAXIT,M1,M2,X0) also returns the
%   iteration number at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = SQMR(A,B,TOL,MAXIT,M1,M2,X0) also returns a
%   vector of the residual norms at each iteration, including NORM(B-A*X0).
%
%   REFERENCES:
%   [1] Freund, Roland W. and Nachtigal, Noel M., "A new Krylov-Subspace Method
%   for Symmetric Indefinite Linear Systems", Proceedings of the 14th IMACS World
%   Congress on Computational and Applied Mathematics, Ed. W. F. Ames, 
%   pp. 1253-1256, IMACS, 1994.
%   [2] Freund, Roland W., "Preconditioning of symmetric but highly indefinite 
%   linear systems", 15th IMACS World Congress on Scientific Computation Modelling
%   and Applied Mathematics, Vol. 2 Numerical Mathematics, Ed. A. Sydow, pp. 551-556,
%   Wissenschaft und Technik Verlag, 1997.

%   Han Chen, 2004.
%   $ Last modified: 2004/03/07 17:03:14 $

% Check for an acceptable number of input arguments
if nargin < 2
   error('Not enough input arguments.');
end

% Determine whether A is a matrix, a string expression,
% the name of a function or an inline object.
[atype,afun,afcnstr] = iterchk(A);
if isequal(atype,'matrix')
   % Check matrix and right hand side vector inputs have appropriate sizes
   [m,n] = size(A);
   if (m ~= n)
      error('Matrix must be square.');
   end
   if ~isequal(size(b),[m,1])
      es = sprintf(['Right hand side must be a column vector of' ...
            ' length %d to match the coefficient matrix.'],m);
      error(es);
   end
else
   m = size(b,1);
   n = m;
   if (size(b,2) ~= 1)
      error('Right hand side must be a column vector.');
   end
end

% Assign default values to unspecified parameters
if nargin < 3 | isempty(tol)
   tol = 1e-6;
end
if nargin < 4 | isempty(maxit)
   maxit = min(n,20);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                      % Norm of rhs vector, b
if (n2b == 0)                       % if    rhs vector is all zeros
   x = zeros(n,1);                  % then  solution is all zeros
   flag = 0;                        % a valid solution has been obtained
   relres = 0;                      % the relative residual is actually 0/0
   iter = 0;                        % no iterations need be performed
   resvec = 0;                      % resvec(1) = norm(b-A*x) = norm(0)      
   if (nargout < 2)
      itermsg('sqmr',tol,maxit,0,flag,iter,NaN);
   end
   return
end

if ((nargin >= 5) & ~isempty(M1))
   existM1 = 1;
   [m1type,m1fun,m1fcnstr] = iterchk(M1);
   if isequal(m1type,'matrix')
      if ~isequal(size(M1),[m,m])
         es = sprintf(['Preconditioner must be a square matrix' ...
               ' of size %d to match the problem size.'],m);
         error(es);
      end      
   end   
else
   existM1 = 0;
   m1type = 'matrix';
end

if ((nargin >= 6) & ~isempty(M2))
   existM2 = 1;
   [m2type,m2fun,m2fcnstr] = iterchk(M2);
   if isequal(m2type,'matrix')
      if ~isequal(size(M2),[m,m])
         es = sprintf(['Preconditioner must be a square matrix' ...
               ' of size %d to match the problem size.'],m);
         error(es);
      end
   end
else
   existM2 = 0;
   m2type = 'matrix';
end

if ((nargin >= 7) & ~isempty(x0))
   if ~isequal(size(x0),[n,1])
      es = sprintf(['Initial guess must be a column vector of' ...
            ' length %d to match the problem size.'],n);
      error(es);
   else
      x = x0;
   end
else
   x = zeros(n,1);
end

if ((nargin > 7) & isequal(atype,'matrix') & ...
      isequal(m1type,'matrix') & isequal(m2type,'matrix'))
   error('Too many input arguments.');
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
if isequal(atype,'matrix')
   r = b - A * x;                  % Zero-th residual
else
   r = b - iterapp(afun,atype,afcnstr,x,varargin{:});
end
normr = norm(r);                   % Norm of residual

if (normr <= tolb)                 % Initial guess is a good enough solution
   flag = 0;
   relres = normr / n2b;
   iter = 0;
   resvec = normr;
   if (nargout < 2)
      itermsg('sqmr',tol,maxit,0,flag,iter,relres);
   end
   return
end

resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin

if existM1
   if isequal(m1type,'matrix')
      t = M1 \ r;
   else
      t = iterapp(m1fun,m1type,m1fcnstr,r,varargin{:});
   end
   if isinf(norm(t,inf))
      flag = 2;
      relres = normr/n2b;
      iter = 0;
      resvec = normr;
      if nargout < 2
         itermsg('sqmr',tol,maxit,0,flag,iter,relres);
      end
      return
   end
else
   t = r;
end
tau = norm(t);

if existM2
   if isequal(m2type,'matrix')
      q = M2 \ t;
   else
      q = iterapp(m2fun,m2type,m2fcnstr,t,varargin{:});
   end
   if isinf(norm(q,inf))
      flag = 2;
      relres = normr/n2b;
      iter = 0;
      resvec = normr;
      if nargout < 2
         itermsg('sqmr',tol,maxit,0,flag,iter,relres);
      end
      return
   end
else
   q = t;
end

thet = 0;
rho = r' * q;

stag = 0;                          % stagnation of the method

% loop over maxit iterations (unless convergence or failure)

for i = 1 : maxit
   if isequal(atype,'matrix')
      t = A * q;
   else
      t = iterapp(afun,atype,afcnstr,q,varargin{:});         
   end
   sigma = q' * t;
   if sigma == 0 | isinf(sigma)
      flag = 4;
      break
   end

   alpha = rho/sigma;
   r = r - alpha*t;
   
   if existM1
      if isequal(m1type,'matrix')
         t = M1 \ r;
      else
         t = iterapp(m1fun,m1type,m1fcnstr,r,varargin{:});
      end
      if isinf(norm(t,inf))
         flag = 2;
         break
      end
   else
      t = r;
   end
   if tau == 0 | isinf(tau)
       flag = 4;
       break;
   end
   thet1 = thet;
   thet = norm(t) / tau;
   c = 1 / sqrt(1 + thet^2);
   if c == 0
       flag = 4;
       break;
   end
   tau = tau * thet * c;
   if i == 1
       d = c^2 * alpha * q;
   else
       d = c^2 * thet1^2 * d + c^2 * alpha * q;
   end
   
   % Check for stagnation of the method
   stagtest = zeros(n,1);
   ind = (x ~= 0);
   stagtest(ind) = d(ind) ./ x(ind);
   stagtest(~ind & d ~= 0) = Inf;
   if norm(stagtest,inf) < eps
      stag = 1;
   end
   
   x = x + d;                       % form the new iterate
   
   if isequal(atype,'matrix')
       normr = norm(b - A * x);
   else
       normr = norm(b - iterapp(afun,atype,afcnstr,x,varargin{:}));
   end
   resvec(i+1) = normr;
   
   if normr <= tolb                 % check for convergence
      flag = 0;
      iter = i;
      break
   end
   
   if stag == 1
      flag = 3;
      break
   end
   
   if normr < normrmin              % update minimal norm quantities
      normrmin = normr;
      xmin = x;
      imin = i;
   end
   
   if rho == 0
       flag = 4;
       break;
   else
       rho1 = rho;
       if existM2
           if isequal(m2type,'matrix')
               u = M2 \ t;
           else
               u = iterapp(m2fun,m2type,m2fcnstr,t,varargin{:});
           end
           if isinf(norm(u,inf))
               flag = 2;
               break
           end
       else
           u = t;
       end
       rho = r' * u;
       beta = rho / rho1;
       q = u + beta * q;
   end
       
end                                % for i = 1 : maxit

% returned solution is first with minimal residual
if flag == 0
   relres = normr / n2b;
else
   x = xmin;
   iter = imin;
   relres = normrmin / n2b;
end

% truncate the zeros from resvec
if flag <= 1 | flag == 3
   resvec = resvec(1:i+1);
else
   resvec = resvec(1:i);
end

% only display a message if the output flag is not used
if nargout < 2
   itermsg('sqmr',tol,maxit,i,flag,iter,relres);
end
