function [ q,ea,iter ] = romberg( f,a,b,es,maxit,varargin )
% simp38: composite simpson 3/8 rule integrator [ (f+3)3n ]
%  [ q,ea,iter ] = romberg( f,a,b,es,maxit,varargin )
%       uses the Romberg iterative integration technique:
%       trapezoidal approximations with increasing segments,
%       followed by repeated Richardson extrapolation using past results
% requires:
%   function I = trap(f,a,b,n,varargin);
% input:
%   f = function to numerically integrate
%   a,b = integration limits
%   es = desired relative error
    esDEF = 0.000001;   % default relative error
%   maxit = maximum allowable iterations
    maxitDEF = 50;  % default maximum number of iterations
%   p1,p2,... = any additional parameters used by f
% output:
%   q = integral estimate
%   ea = approximate relative error (%)
%   iter = number of iterations
%##########################################################################
% Pseudo Code: (flops per section)
%   Total Flops: 
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations: 
%   ====
%   Main Algorithm: 
%   ####
%##########################################################################
% Input Format Check:

if nargin<3,error('Need function and integration limits');end
if nargin<4 || isempty(es), es=esDEF;end
if nargin<5 || isempty(maxit), maxit=maxitDEF;end

%==========================================================================
% Variable Declarations:
%
% Total Flops: trap(f,n)

n = 1;  % number of integration segments for each step, starting with one
I(1,1) = trap(f,a,b,n,varargin{:}); % grid of guesses beginning with trap
iter = 0;   % iteration counter

%==========================================================================
% Main Algorithm:
%
% Total Flops:

% while within our maximum number of iterations,
% approximate the integral using an increasing number of segments
while iter<maxit
    iter = iter +1; % incremement the iteration counter
    n = 2^iter; % double the number of segments in the approximation
    
    % Compute and add the new approximation
    I(iter+1,1) = trap(f,a,b,n,varargin{:});
    
    % use the new approximation to compute a new set of refined ones
    % using successive Richardson extrapolations
    for k = 2:(iter+1)
        j = 2+iter-k;
        I(j,k) = (4^(k-1)*I(j+1,k-1)-I(j,k-1))/(4^(k-1)-1);
    end
    
    % compute the relative error and break if we are within tolerance
    ea = abs((I(1,iter+1)-I(2,iter))/I(1,iter+1))*100;
    if ea<=es, break;end
end
q = I(1,iter+1);    % return the best approximation

%##########################################################################