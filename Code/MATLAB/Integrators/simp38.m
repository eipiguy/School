function I = simp38( f,a,b,n,varargin )
% simp38: composite simpson 3/8 rule integrator [ (f+3)3n ]
%  I = simp38( f,a,b,n,varargin )
%       estimates the integral of f from a to b
%       using simposins 3/8 rule with n segments 
%       (picking 4 points from each with matching endpoints)
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
%   Total Flops: (f+3)3n
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations: (f+3)
%   ====
%   Main Algorithm: (f+3)3n -(f+3)
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
% Total Flops: (1 +2 +f) = (f+3)

x = a;                  % x value "buffer" for each segment
n = 3*n;                % double the number of segments for the algorithm
h = (b-a)/n;            % segment width
s = f(a,varargin{:});   % sum of f values, added in segments and scaled
% 1 flops for doubling the number of segments
% 2 flops to calculate the segment width
% f flops for calculating f at the first point

j = 2;          % current factor in use
c = [2;3;3];    % constants for weighting the algorithm correctly
l = length(c);  % length of coefficients vector

%==========================================================================
% Main Algorithm:
%
% Total Flops: sum_{i=2}^{3n}[ 1 +f +2 ] +f +1 +2 = (3n-1)(f+3)
%   = (f+3)3n -(f+3)

for i=2:n
    x = x + h;
    s = s + (c(j)*f(x,varargin{:}));
    switch j
        case l
            j = 1;
        otherwise
            j = j+1;
    end
end
% (n-2) cycles in i
%   1 flop each for shifting the x value to the next segment boundary
%   f flops for computing the new f value,
%   2 for scaling by the correct factor and adding to the partial sum

s = s + f(b,varargin{:});
I = s*h*3/8;
% (f +1) flops for computing the last value of f and adding it to the sum
% 2 flops for scaling the final sum to the proper value

%##########################################################################