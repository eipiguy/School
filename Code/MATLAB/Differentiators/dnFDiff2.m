function df = dnFDiff2( n,f,x,h,varargin )
% dnFDiff: 
%   df = dnFDiff( n,f,x,h,varargin )
% input:
%   n = the order of the derivative approximation
%   f = function to differentiate
%   x = x value to evaluate the derivative at
%   h = step size
%   p1,p2,... = additional parameters used by dydt
% output:
%   df = approximate derivative of f at x

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

if nargin<4 || isempty(h),h=1e-6;end
if nargin<3,error('Need d-order, function, and evaluation point');end

%==========================================================================
% Variable Declarations:

xspan = x:h:(x+(n+1)*h);
fspan = f(xspan);

c = zeros(n+2,1);
for i=0:n
    c(i+1) = ((-1)^(n+i))*((factorial(n+1)*nchoosek(n,i)) +nchoosek(n+1,i));
end
c(n+2) = ((-1)^(n+1+i));
disp(xspan);
disp(fspan);
disp(c);

%==========================================================================
% Main Algorithm:

df = (fspan*(c))./( factorial(n+1) *(h.^(n)) );

%##########################################################################
