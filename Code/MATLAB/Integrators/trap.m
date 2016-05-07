function I = trap( f,a,b,n,varargin )
% trap: composite trapezoidal rule integrator [(f+3)n+(f+4)]
%  I = trap( f,a,b,n,varargin )
%       estimates the integral of f from a to b
%       using the trapezoidal rule with n segments
% input:
%   f = function to numerically integrate
%   a,b = integration limits
%   n = number of segments in integral approximation
    nDEF = 100; % default number of segments
%   p1,p2,... = any additional parameters used by f
% output:
%   I = integration estimate
%##########################################################################
% Pseudo Code: (flops per section)
%   Total Flops: (f+3)n+(f+4)
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations: (f+2)
%   ====
%   Main Algorithm: (f+3)n
%   ####
%##########################################################################
% Input Format Check:

% make sure we have a function to integrate and limits of the integral
if nargin < 3, error('Need a function and limits of the integral.');end

% make sure the larger number in the integral is on the top
if ~(b>a), error('Upper bound must be greater than lower bound.');end

% set the number of segments to the default
% if they weren't already input
if nargin<4 || isempty(n),n=nDEF;end

%==========================================================================
% Variable Declarations:
%
% Total Flops: ( f +2 )

x = a;                  % x value "buffer" for each segment
h = (b-x)/n;            % segment width
s = f(x,varargin{:});   % sum of f values, added in segments and scaled
% 2 flops to calculate the segment width
% f flops for calculating f at the first point

%==========================================================================
% Main Algorithm:
%
% Total Flops: sum_{i=1}^{n-1}[ 1 +f +2 ] + ( f +1 ) +2
%   = (f+3)n

for i=1:(n-1)
    x = x + h;
    s = s + 2*f(x,varargin{:});
end
% (n-1) cycles in i
%   1 flop each for shifting the x value to the next segment boundary
%   f flops for computing the new f value,
%   2 for doubling it, and then adding it to the partial sum

s = s + f(b,varargin{:});
I = s*h/2;
% (f +1) flops for computing the last value of f and adding it to the sum
% 2 flops for scaling the final sum to the proper value

%##########################################################################