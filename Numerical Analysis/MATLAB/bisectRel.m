function [root,fx,ea,iter] = bisectRel(func,xl,xu,es,maxit,varargin)
% bisect: root locations
%   [root,fx,ea,iter] = bisect(func,x1,xu,es,maxit,p1,p2,...)
%       uses bisection to find a root of func 
%       within an relative error margin (% of estimate)
% input:
%   func = name of function
%   xl,xu = lower and upper guesses
%   es = desired absolute error (default = 0.0001)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = additional parameters used by func
% output:
%   root: real root
%   fx = funtion value at root
%   ea = approximate relative error
%   iter = number of iterations

% make sure a function and upper/lower guesses were given
if nargin < 3,error('at least 3 input arguments required'),end

% make sure there is an odd number of roots within the guesses
% (hopefully only one)
test = func(xl,varargin{:})*func(xu,varargin{:});
if test > 0,error('no sign change'),end

% set default desired error and maximum iterations
if nargin < 4 | isempty(es),es=0.0001;end
if nargin < 5 | isempty(maxit),maxit=50;end

% initialize iterations, root guess, and resulting error
iter = 0; xr = xl; ea = 100;
% while we have not reached the maximum iterations needed
while(1)
    % save the last root guess, and pick the new one, and step
    xrold = xr;
    xr = (xl + xu)/2;
    iter = iter + 1;
    
    % if we are not at an exact root,
    % compute the current relative error
    if xr ~= 0,ea = abs((xr - xrold)/xr)*100;end
    
    % check to see which half of the current interval
    % has an odd number of roots (hopefully just one root)
    % and set it as the new interval for the next iteration
    test = func(xl,varargin{:})*func(xr,varargin{:});
    if test < 0
        xu = xr;
    elseif test > 0
        xl = xr;
    else
        ea = 0;
    end
    
    % if we've hit the maximum number of iterations,
    % stop the loop and set the values to return
    if ea <= es | iter >= maxit,break,end
end
root = xr; fx = func(xr,varargin{:});