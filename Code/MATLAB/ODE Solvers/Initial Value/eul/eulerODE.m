function [ t,y ] = eulerODE( dydt,tspan,y0,h,varargin )
% eulerODE: solves a first order ODE using Euler's method
%   [ t,y ] = eulerODE( dydt,tspan,y0,h,varargin )
%       Uses Euler's method to numerically integrade a first order ODE
% input:
%   dydt = function that evaluates the derivative
%   tspan = [ti,tf]
%       ti = initial time
%       tf = final time
%   y0 = initial value of dependent variable
%   h = step size
%   p1,p2,... = additional parameters used by dydt
% output:
%   t = vector of independent variable
%   y = vector of solution for dependent variable
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

if nargin<4,error('At least 4 arguments required.');end

ti = tspan(1);tf = tspan(2);
if ~(tf>ti),error('Upper limit must be greater than the lower'),end;

%==========================================================================
% Variable Declarations:

t = (ti:h:tf)'; n = length(t);
% if necessary, add an additional value of t
% so the range goes from t = ti to tf
if t(n)<tf
    t(n+1) = tf;
    n = n+1;
end

%==========================================================================
% Main Algorithm:

y = ones(n,1)*(y0);    % preallocate y to improve efficiency
for i = 1:(n-1)         % Step through Euler's method
    y(i+1,:) = y(i,:) + ((t(i+1)-t(i))*dydt(t(i),y(i,:),varargin{:}))';
end
%==========================================================================
%% Plot Results

m = length(y0);
sp = 1;     % Subplot index

cla
hold off
for i=1:m
    subplot(m,m,sp),plot(t,y(:,i));
    xlabel('t');
    ylabel(sprintf('y(%d)',i));
    sp = sp+1;
end

for i=1:m-1
    for j=i+1:m
        subplot(m,m,sp),plot(y(:,i),y(:,j));
        xlabel(sprintf('y(%d)',i));
        ylabel(sprintf('y(%d)',j));
        sp = sp+1;
    end
end
%##########################################################################
end