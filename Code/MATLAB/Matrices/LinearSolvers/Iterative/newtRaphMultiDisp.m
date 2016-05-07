function [ x,f,ea,iter ] = newtRaphMultiDisp( func,x0,es,maxit,varargin )
% newtRaphMultiDisp: 
%   [ x,f,ea,iter ] = newtRaphMultiDisp( func,x0,es,maxit,p1,p2,... )
%       Uses the multivariable Newton-Raphson method
%       to approximate the roots of a nonlinear, multivariable function
%       while displaying results at each step
% input:
%   func = name of an input multivariable function that returns
%       a functional value and the Jacobian
%   x0 = initial guess for root
%   es = desired percent relative error (default = 0.0001%)
%   maxit = maximum allowed iterations (default = 50)
%   p1,p2,... = additional parameters used by func
% output:
%   x = solution root vector
%   f = func evaluated at the given root vector
%   ea = approximate percent relative error
%   iter = number of iterations
%##########################################################################
% Pseudo Code:
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Main Algorithm:
%   ####
%##########################################################################
% Input Format Check:

% Check to make sure there is a function and a rhs vector.
if nargin<2, error('function and initial guess required');end

% Set defaults for maximum iterations and relative error
% if they aren't explicitly given
if nargin<4 || isempty(maxit), maxit=50; end
if nargin<3 || isempty(es), es=0.0001; end

%==========================================================================
% Variable Declarations:

% initialize the iteration counter and the initial guess
iter = 0;
x = x0;
display(x);

%==========================================================================
% Main Algorithm:

while(1)
    % The function should return its own Jacobian at each point
    [J,f] = func(x,varargin{:});
    
    % Solve for J^{-1}(f) and iterate using Newt-Raph algorithm
    dx = J\f;
    display([f,J]);
    display([x,dx]);
    
    x = x -dx;
    display(x);
    
    % Increment the iteration counter
    iter = iter+1;
    
    % Break if we are within our tolerances
    ea = max(abs(dx./x))*100;
    display([iter,ea]);
    
    if iter >= maxit || ea <= es, break, end
end

%##########################################################################
end

