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

m = length(y0);                 % Number of variables in the system

ti = tspan(1);tf = tspan(2);    % End points for independent variable
t = (ti:h:tf)';                 % Span of independents variable
n = length(t);                  % Number of steps

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

% Clear all current figures and turn off hold
cla
hold off

% Set up the grid of sub plots for each of the time and phase plots.
for i=1:(m^(2))
    sub(i) = subplot(m,m,i);
    set(sub(i),'Visible','off');
end

% Place each time plot solution in the left column.
for i=1:m
    subplot(sub(1+(i-1)*m)),plot(t,y(:,i));
    xlabel('t');
    ylabel(sprintf('y(%d)',i));
end

% Place phase graphs next to their corresponding time solutions
for i=1:m-1
    for j=i+1:m
        subplot(sub(((i-1)*m)+j)),plot(y(:,j),y(:,i));
        xlabel(sprintf('y(%d)',j));
        ylabel(sprintf('y(%d)',i));
    end
end
%##########################################################################
end