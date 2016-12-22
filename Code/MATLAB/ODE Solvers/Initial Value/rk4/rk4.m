function [ tOut,yOut ] = rk4( dydt,tspan,y0,h,varargin )
%% rk4ODEsysVid: numerically solves a system of first order ODEs using RK4
%   [ t,y ] = rk4sys( dydt,tspan,y0,h,varargin )
%       Uses the 4th order Runga Kutta method
%       to numerically solve a system of first order ODEs
%       and output solutions in row vectors stacked by time in a column,
%       given a row vector with the initial positions of each component.
%##########################################################################
%% Input:
%   dydt = differential equation, first two inputs must be as ( t, y, ... )
%           t = independent column vector (t)
%           y = dependent row vector (y)
%               = [ y(1), y(2), ... ]
%           varargin = any additional parameters needed
%       output derivative values must be row vectors
%           dydt = [ dydt(1), dydt(2), ... ]
%   tspan = [ti,tf]
%       ti = initial time
%       tf = final time
%   OR tspan = [ti,t1,t2...,tf] 
%           where each ti is a time needed in the approximation
%   y0 = row vector for initial values of dependent variables
%       = [ y0(1), y0(2), ... ]
%   h = time interval step size
%   varargin = p1,p2,... = additional parameters used by dydt
%==========================================================================
%% Output:
%   tOut = row vector of times for each column of the solution
%           = [ tOut(1), tOut(2), ... ]
%   yOut = row vector with columns for each solution at a given time
%           = [ yOut(:,1), yOut(:,2), ... ]
%##########################################################################
%% Pseudocode:
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Main Algorithm:
%   ####
%##########################################################################
%% Input Format Check:

% Make sure all inputs are given,
if nargin<4
    error('ODE, time span, initial conditions, and step size required.');
end
% Make sure initial and final times are in increasing order
if any(diff(tspan)<=0),error('tspan not in ascending order');end

%==========================================================================
%% Variable Declarations:

m = length(y0);                 % Number of variables in the system
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
    t = tspan; % Otherwise the times to approximate are given explicitly
end

% Set up the other initial parameters for the method
tt = ti; y(1,:) = y0;
np = 1; tOut(np) = tt; yOut(np,:) = y(1,:);
i = 1;

%==========================================================================
%% Interval Algorithm:

% For each given time point:
while(1)
    tend = t(np+1);         % Set the current sub-interval's end-point
    hh = t(np+1) - t(np);   % Match the current sub-intervals's step size
        % 1 flop for setting the sub-step size
    
    % If the time intervals were given explicitly, 
    % the sub-interval may have a larger step size than the input step.
    % If so, we adjust this sub-step size to the smaller value.
    if hh>h, hh=h; end
    
    % Approximate the function value for the end of the sub-interval:
    while(1)
        % If the sub-step size still overshoots the needed sub-interval,
        % then we shrink it to fit again.
        if tt+hh>tend, hh = tend-tt;end
            % 1 or 2 flops, 
            % one for checking the interval and one for adjusting it
        
        % Find the slope at the beginning of the sub-interval (row vector)
        k1 = dydt(tt,y(i,:),varargin{:});
            % dydt flops for computing the slope
        
        % Project k1 linearly to get an estimate for the function
        % at the midpoint of the sub-interval, (row vector)
        ymid = y(i,:) + (k1.*hh./2);        
        % and calculate the slope at this midpoint. (another row vector)
        k2 = dydt(tt+(hh/2),ymid,varargin{:});
            % 3m flops for projecting
            % 3 flops for setting the time value
            % dydt flops for computing the slope
        
        % Using the slope at the midpoint k2, linearly
        % recalculate the projection to the midpoint
        ymid = y(i,:) + (k2.*hh./2);        
        % and find another slope approximation k3 for the midpoint.
        k3 = dydt(tt+hh/2,ymid,varargin{:});
            % 3m flops for projecting
            % 3 flops for setting the time value
            % dydt flops for computing the slope
        
        % Using k3, linearly re-project all the way 
        % to the end of the interval
        yend = y(i,:) + (k3.*hh);
        % and approximate the slope there.
        k4 = dydt(tt+hh,yend,varargin{:});
            % 3m flops for projecting
            % 2 flops for setting the time value
            % dydt flops for computing the slope
        
        % Use the weighted RK4 formula to refine the approximation.
        phi = ( k1 + 2*(k2+k3) + k4 )/6;
            % 5m flops for the final slope approximation
        
        % With the rk4 slope approximation, project linearly one more time 
        % to approximate the function at the end of the sub-interval,
        % and take this as the value for start of the next sub-interval
        y(i+1,:) = y(i,:) + (phi*hh);
            % 2m flops for projecting to the next step
        
        % Move to the next sub-interval,
        % and break once we reach the next given time point
        tt = tt+hh;
        i = i+1;
        if tt>=tend,break,end
            % 1 flop for moving to the next step
    end
    
    % Once we have reached the end of each time interval,
    % we record the values for output, and move to the next one.
    np = np+1; tOut(np) = tt; yOut(np,:) = y(i,:);
    
    % Once we hit the end of the parent interval, we break the loop.
    if tt>=tf,break,end
end

%##########################################################################
end