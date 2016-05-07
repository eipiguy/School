function [ t,y ] = eulExODEsysD( dydt,tspan,y0,h,varargin )
%% eulExODEsysD: solves a system first order ODEs using Euler's method
%   [ t,y ] = eulExODEsysD( dydt,tspan,y0,h,varargin )
%       Uses the explicit Euler's method to 
%       approximate a system first order ODEs
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
%% Pseudo Code:
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Main Algorithm:
%   ====
%   Plot Results:
%   ####
%##########################################################################
%% Input Format Check:

% Make sure all inputs are given.
if nargin<4
    error('ODE, time span, initial conditions, and step size required.');
end

% Make sure initial and final times are in increasing order
if any(diff(tspan)<=0),error('tspan not in ascending order');end

%==========================================================================
%% Variable Declarations:

ti = tspan(1);tf = tspan(2);
t = (ti:h:tf)'; n = length(t);

% if necessary, add an additional value of t
% so the range goes from t = ti to tf
if t(n)<tf
    t(n+1) = tf;
    n = n+1;
end

%==========================================================================
%% Main Algorithm:

y = ones(n,1)*(y0);    % preallocate y to improve efficiency
for i = 1:(n-1)         % Step through Euler's method
    y(i+1,:) = y(i,:) + ((t(i+1)-t(i))*dydt(t(i),y(i,:),varargin{:}))';
end

%==========================================================================
%% Plot Results

m = length(y0); % Numer of independent variables
sp = 1;         % Subplot index

cla
hold off
for i=1:m
    subplot(m-1,m,sp),plot(t,y(:,i));
    xlabel('t');
    ylabel(sprintf('y(%d)',i));
    sp = sp+1;
end

for i=1:m-1
    for j=i+1:m
        subplot(m-1,m,sp),plot(y(:,i),y(:,j));
        xlabel(sprintf('y(%d)',i));
        ylabel(sprintf('y(%d)',j));
        sp = sp+1;
    end
end
%##########################################################################
end