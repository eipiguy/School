function df = dnBDiff( n,f,x,h,varargin )
% dnBDiff: 
%   df = dnBDiff( n,f,x,h,varargin )
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

if nargin<4 || isempty(h),h=.5e-6;end
if nargin<3,error('Need d-order, function, and evaluation point');end

%==========================================================================
% Variable Declarations:

xspan = x:-h:(x-n*h);
fspan = f(xspan);

c = zeros(n+1,1);
for i=0:n
    c(i+1) = ((-1)^(i))*nchoosek(n,i);
end
disp(xspan);
disp(fspan);
disp(c);

%==========================================================================
% Main Algorithm:

df = (fspan*c)./(h.^(n));

%##########################################################################