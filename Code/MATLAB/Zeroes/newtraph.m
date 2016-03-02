function [root,ea,iter] = newtraph(func,dfunc,xr,es,maxit,varargin)
% newtraph: Newton-Raphson root locations
%   [root,ea,iter] = newtraph(func,dfunc,xr,es,maxit,p1,p2,...):
%       uses Newton-Raphson method to find the root of func
% input:
%   func = name of function
%   dfun = name of derivative of function
%   xr = initial guess
%   es = desired relative error (default = 0.0001%)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = additional parameters used by function
% output:
%   root: real root
%   ea = approximate relative error (%)
%   iter = number of iterations

% make sure a function, its derivative, and an initial guess were given
if nargin < 3,error('at least 3 input arguments required'),end

% set default desired error and maximum iterations
if nargin < 4 || isempty(es),es=0.0001;end
if nargin < 5 || isempty(maxit),maxit=50;end

iter = 0;   % initialize iterations

% while we have not reached the maximum iterations needed
while(1)
    xrold = xr;                     % save the last root guess,
    xr = xr - (func(xr)/dfunc(xr)); % pick the new one,
    iter = iter + 1;                % and step the iteration counter
    
    % if we are not at an exact root,
    % compute the current relative error
    if xr ~= 0,ea = abs((xr - xrold)/xr)*100;end
    
    % if we've reached the solution within desired error
    % or hit the maximum number of iterations,
    % stop the loop and set the values to return
    if ea <= es || iter >= maxit,break,end
end
root = xr;