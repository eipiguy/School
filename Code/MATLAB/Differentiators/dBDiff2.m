function df = dBDiff2( n,f,x,h,varargin )
% dBDiff2: 
%   df = dBDiff2( n,f,x,h,varargin )
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

xspan = x:-h:(x-((n+1)*h));
fspan = f(xspan);

switch n
    case 1
        c = [3,-4,1]./2;
    case 2
        c = [2,-5,4,-1];
    case 3
        c = [5,-18,24,-14,3]./2;
    case 4
        c = [3,-14,26,-24,11,-2];
    otherwise
        error('Can only compute up to 4th derivative.');
end
disp(xspan);
disp(fspan);
disp(c);

%==========================================================================
% Main Algorithm:

df = (fspan*(c'))./(h.^(n));

%##########################################################################
