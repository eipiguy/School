function [ tp,yp ] = rk4ODEsysDisp( dydt,tspan,y0,h,varargin )
%% rk4ODEsysDisp: solves a system of first order ODEs, and plots the result
%   [ t,y ] = rk4sys( dydt,tspan,y0,h,varargin )
%       Uses the 4th order Runga Kutta method
%       to numerically integrade a system of first order ODEs,
%       and then plot the results
% input:
%   dydt = differential equation with independent(t) & dependent(y) inputs
%   tspan = [ti,tf]
%       ti = initial time
%       tf = final time
%   OR tspan = [ti,t1,t2...,tf] points to approximate the function at
%   y0 = initial value of dependent variable
%   h = step size
%   p1,p2,... = additional parameters used by dydt
% output:
%   tp = vector of independent variable
%   yp = vector of solution for dependent variable

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

% Make sure all inputs are given,
if nargin<4,error('At least 4 arguments required.');end

% Make sure initial and final times are in increasing order
if any(diff(tspan)<=0),error('tspan not in ascending order');end

%==========================================================================
%% Variable Declarations:

n = length(tspan);              % Number of steps between endpoints
ti = tspan(1); tf = tspan(n);   % Set variables first time and last time

% Setting up the steps for the time span:
% If the time span only contains an initial and final time,
if n==2
    t = (ti:h:tf)'; n = length(t);  % fill in the steps between.
    
    % If the last element didn't reach the final time
    if t(n)<tf
        t(n+1) = tf;    % Add an extra step
        n = n+1;
    end
else
    t = tspan; % Otherwise the times to approximate at are given explicitly
end

% Set up the other initial parameters for the method
tt = ti; y(1,:) = y0;
np = 1; tp(np) = tt; yp(np,:) = y(1,:);
i = 1;

%==========================================================================
%% Main Algorithm:

% For each integral step in the approximation:
while(1)
    tend = t(np+1);         % Set the current interval's end-point
    hh = t(np+1) - t(np);   % and the current interval's step size
    
    % If the time steps were given explicitly, 
    % this may be larger than the given step size.
    % If so, we shrink it down and compute intermediate values.
    if hh>h, hh= h;end
    
    % Aprroximate the function value for the start of the next step:
    while(1)
        % If the current step size overshoots the current interval,
        % then we also chop it down to fit
        if tt+hh>tend, hh = tend-tt;end
        
        % Start with the slope at the beginning of the current interval
        k1 = dydt(tt,y(i,:),varargin{:})';
        
        % Project k1 to get an estimate for the value of the function
        % at the midpoint of the current interval,
        ymid = y(i,:) + k1.*hh./2;        
        % and calculate the slope at this midpoint.
        k2 = dydt(tt+hh/2,ymid,varargin{:})';
        
        % Using the new slope k2, 
        % recalculate the projection at the midpoint
        ymid = y(i,:) + k2.*hh./2;        
        % and find another new slope k3.
        k3 = dydt(tt+hh/2,ymid,varargin{:})';
        
        % Using the last slope k3,
        % Re-project all the way to the end of the interval
        yend = y(i,:) + k3.*hh;
        % and approximate the slope here.
        k4 = dydt(tt+hh,yend,varargin{:})';
        
        % Use the weighted RK4 formula to refine the approximation.
        phi = (k1 +2*(k2+k3) +k4)/6;
        
        % With the final slope approximation, project one more time 
        % to find the value of the function for the start of the next step
        y(i+1,:) = y(i,:) +phi*hh;
        
        % Move to the next intermediate step in the approximation,
        % and break once we reach the next given time to check the function
        tt = tt+hh;
        i = i+1;
        if tt>=tend,break,end
    end
    
    % Once we have reached a specified time in tspan,
    % we record the values for output, and move to the next one.
    np = np+1; tp(np) = tt; yp(np,:) = y(i,:);
    
    % Once we hit the end of the full interval, we break the loop.
    if tt>=tf,break,end
end
%##########################################################################
end