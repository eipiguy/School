function [ p ] = hornetPolynomial( a,x )
% hornetPolynomial: value of a polynomial
%   [ f ] = hornetPolynomial( a,x )
%       uses the hornet alogrithm to compute the value of a polynomial
%       with the given coefficients a at the point x
% input:
%   a = vector containing the coefficients of the polynomial
%       (ordered by increasing power)
%   x = point at which to evaluate the given polynomial
% output:
%   f = functional value of the polynomial at the given point
%
% This algorithm follows the idea of factoring a polynomial
% to give a recursive formula that avoids redundant computations of x^n
% For example a cubic can be written as:
%   p(x) = a1 + a2*x + a3*x^2 + a4*x^3 = a1 + x*(a2 + x*(a3 + x*a4))

% Initialize the function by setting it to the leading power's constant 
p = a(length(a));

% Multiply the leading term by x and add the constant for the 
% next smallest power until all the constants have been exhausted
for i= length(a)-1:-1:1
    p = p*x + a(i);
end

end

