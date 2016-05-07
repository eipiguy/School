function [ p ] = hornet( x,a )
% hornet: recursive polynomial evaluation [ 2n -2 ]
%   [ f ] = hornet( x,a )
%       uses the hornet alogrithm to compute the value of a polynomial
%       with the given coefficients a (ordered by increasing power)
%       at the point x
% input:
%   x = point at which to evaluate the polynomial
%   a = coefficients of the polynomial (ordered by increasing power)
% output:
%   f = functional value of the polynomial at the given point
%
% This algorithm follows the idea of factoring a polynomial
% to give a recursive formula that avoids redundant computations of x^n
% For example a cubic can be written as:
%   p(x) = a1 + a2*x + a3*x^2 + a4*x^3 = a1 + x*(a2 + x*(a3 + x*a4))
%##########################################################################
% Pseudo Code: (flops per section)
%   Total Flops: 2n -2
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Main Algorithm: 2n -2
%   ####
%##########################################################################
% Input Format Check:

%==========================================================================
% Variable Declarations:

n = length(a); % largest power of the polynomial

% Initialize the function by setting it to the leading power's constant 
p = a(n);   % recursively scaled polynomial value

%==========================================================================
% Main Algorithm:
%
% Total Flops: sum_{i=1}^{n-1}[ 2 ] = 2(n-1) = 2n -2

% Multiply the leading term by x and add the constant for the 
% next smallest power until all the constants have been exhausted
for i= n-1:-1:1
    p = p.*x + a(i);
end
% (n-1) cycles in i
%   2 flops each for scaling the last value and adding in the new constant

end

