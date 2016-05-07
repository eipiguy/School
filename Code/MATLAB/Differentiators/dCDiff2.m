function df = dCDiff2( n,f,x,h,varargin )
% dCDiff2: 
%   df = dCDiff2( n,f,x,h,varargin )
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

if mod(n,2)==1,xspan = (x+(((((n+1)/2)+1)*h))):-h:(x-(((((n+1)/2)+1)*h)));
else
    xspan = (x+(((n/2)+1)*h)):-h:(x-(((n/2)+1)*h));
end
fspan = f(xspan);

switch n
    case 1
        c = [-1,8,0,-8,1]./12;
    case 2
        c = [-1,16,-30,16,-1]./12;
    case 3
        c = [-1,8,-13,0,13,-8,1]./8;
    case 4
        c = [-1,12,-39,56,-39,12,-1]./6;
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
