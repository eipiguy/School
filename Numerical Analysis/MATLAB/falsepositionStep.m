function [xr,fx] = falsepositionStep(func,xl,xu,varargin)
% false position step: finding roots
%   [xr,ea] = falsepositionStep(func,xl,xu,varargin)
%       uses false position algorithm to guess the root
% input:
%   func = name of function
%   xl,xu = lower and upper guesses
%   p1,p2,... = additional parameters used by func
% output:
%   xr = next iterative root guess
%   fx = funtion value at xr

% make sure a function and upper/lower guesses were given
if nargin < 3,error('at least 3 input arguments required'),end

% make sure there is an odd number of roots within the guesses
% (hopefully only one)
test = func(xl,varargin{:})*func(xu,varargin{:});
if test > 0,error('no sign change'),end

% guesses the root using the false position algorithm
% and checks the functional value at the guess
xr = xu - func(xu)*(xl - xu)/(func(xl)-func(xu));
fx = func(xr,varargin{:});

end

