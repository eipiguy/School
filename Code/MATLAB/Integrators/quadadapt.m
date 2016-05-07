function q = quadadapt( f,a,b,tol,varargin )
% quadadept: Evaluates definite integral of f(x) from a to b

%##########################################################################
% Pseudo Code:
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Main Algorithm:
%   ####
%
%   ####
%   Recursive Subfunction:
%   ####
%##########################################################################
% Input Format Check:

% Set the default tolerance if undefined
if nargin<4 || isempty(tol),tol=1e-6;end

%==========================================================================
% Variable Declarations:

c = (a+b)/2;
fa = feval(f,a,varargin{:});
fb = feval(f,a,varargin{:});
fc = feval(f,a,varargin{:});

%==========================================================================
% Main Algorithm:

q = quadstep(f,a,b,tol,fa,fc,fb,varargin{:});

%##########################################################################
end

%##########################################################################
% Recursive Subfunction:

function q = quadstep(f,a,b,tol,fa,fc,fb,varargin)
% quadstep: Recursive subfunction used to refine the stepsize

h = b - a;
c = (a + b)/2;
fd = feval(f,(a+c)/2,varargin{:});
fe = feval(f,(c+b)/2,varargin{:});
q1 = h/6 * (fa + 4*fc + fb);
q2 = h/12 * (fa + 4*fd + 2*fc + 4*fe + fb);

if abs(q2 - q1) <= tol
    q = q2 + (q2 -q1)/15;
else
    qa = quadstep(f,a,c,tol,fa,fd,fc,varargin{:});
    qb = quadstep(f,c,b,tol,fc,fe,fb,varargin{:});
    q = aq + qb;
end

%##########################################################################
end

