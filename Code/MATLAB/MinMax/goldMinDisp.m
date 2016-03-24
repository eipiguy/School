function [x,fx,ea,iter] = goldMinDisp(f,xl,xu,es,maxit,varargin)
% goldsecDisp: golden section search for minimum
%   [root,fx,ea,iter] = goldMinDisp(func,x1,xu,es,maxit,p1,p2,...)
%       uses golden section search to find the minimum of f
%       within some relative error
% input:
%   func = name of function
%   xl,xu = lower and upper guesses
%   es = desired absolute error (default = 0.0001)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = additional parameters used by func
% output:
%   x = location of minimum
%   fx = minimum funtion value
%   ea = approximate relative error (%)
%   iter = number of iterations

% make sure a function and upper/lower guesses were given
if nargin < 3,error('at least 3 input arguments required'),end

% set default desired error and maximum iterations
if nargin < 4 || isempty(es),es=0.0001;end
if nargin < 5 || isempty(maxit),maxit=50;end

% initialize the golden ratio
phi = (1+sqrt(5))/2;

% initialize iterations, and resulting error
iter = 0; ea = 100;
% while we have not reached the maximum iterations needed
while(1)
    % set the stepping distance
    d = (phi-1)*(xu-xl);
    
    % find intermediary points
    x1 = xl+d;
    x2 = xu-d;
    
    % display the current bounds, intermediates,
    % and functional values
    display([xl,xu,x1,x2]);
    display([f(x1),f(x2)]);
    
    % find which intermediate value is closest to the minimum
    % and revise the upper and lower guesses
    % including it as one of the appropriate new end points
    if f(x1,varargin{:}) < f(x2,varargin{:})
        xopt = x1;
        xl = x2;
    else
        xopt = x2;
        xu = x1;
    end
    iter = iter +1; %step the numer of iterations
    
    % if we are not at an exact root,
    % compute the current relative error
    if xopt ~= 0,ea = (2-phi)*abs((xu-xl)/xopt)*100;end
    
    % display the approximate error
    display(ea);
        
    % if we've hit the maximum number of iterations,
    % or gotten within our error bounds,
    % stop the loop and set the values to return
    if ea <= es || iter >= maxit,break,end
end
x = xopt; fx = f(xopt,varargin{:});
