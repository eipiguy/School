function [ tOut,yOut,fireT,threshT ] = ...
    rk4_resThresh( dydt,tspan,y0,h,thresholdMask,threshold,holdValue,holdTime,varargin )
%% rk4_resThresh: solves a system of ODEs with a threshold and hold time
%   [ t,y ] = rk4_resThresh( dydt,tspan,y0,h, ...
%                   thresholdMask,threshold,holdValue,holdTime,varargin )
%       Uses the 4th order Runga Kutta method
%       to numerically solve a system of first order ODEs,
%       and then plot the results in both time and phase.
% input:
%   dydt = differential equation with independent(t) & dependent(y) inputs
%   tspan = [ti,tf]
%       ti = initial time
%       tf = final time
%   OR tspan = [ti,t1,t2...,tf] points to approximate the function at
%   y0 = initial values of dependent variables
%   h = step size
%   thresholdMask = row vector of zeros and ones 
%                   where a one indicates to threshold check
%                   the corresponding variable in the system
%   threshold = threshold firing value
%   holdValue = firing hold state
%   holdTime = time to hold after threshold fire
%   p1,p2,... = additional parameters used by dydt
% output:
%   tOut = vector of independent variables
%   yp = vector of solution for dependent variables
%##########################################################################
%% Pseudo Code:
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Main Algorithm:
%       ----
%       Threshold:
%       ----
%   ====
%   Plot Results:
%   ====
%   Make Video Frames:
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
    t = tspan; % Otherwise the times to approximate at are given explicitly
end

% Set up the other initial parameters for the method:

% Inner solver starting conditions
i=1; tt = ti; y(1,:) = y0;

% Threshold condition initializations
resetTime = zeros(m);
threshT = zeros(n,m);
fireT = zeros(n,m);

%--------------------------------------------------------------------------
%% Initial Threshold Check
for j=1:m
    % If we start above a threshold value, 
    % then we change the initial value to match firing behavior
    if (y0(j)>=threshold) && thresholdMask(j);
        y(i,j) = holdValue;
        threshT(i,j) = true;
        fireT(i,j) = true;
        resetTime(j) = tt + holdTime;
    end
end
%--------------------------------------------------------------------------

% Partition starting conditions
np = 1; tOut(np) = tt; yOut(np,:) = y(1,:);

%==========================================================================
%% Main Algorithm:

% For each given independent value:
while(1)
    tend = t(np+1);         % Set the current interval's end-point
    hh = t(np+1) - t(np);   % and the current interval's step size
        % 1 flop for setting the step size
    
    % If the time intervals were given explicitly, 
    % the current interval may be larger than the given step size.
    % If so, we adjust the step size and compute intermediate values.
    if hh>h, hh= h;end
    
    % Approximate the function value for the start of the next step:
    while(1)
        % If the step size overshoots the current interval,
        % then we also chop it down to fit
        if tt+hh>tend, hh = tend-tt;end
            % 1 or 2 flops, 
            % one for checking the interval and one for adjusting as needed
        
        % Start with the slope at the beginning of the step
        k1 = dydt(tt,y(i,:),varargin{:});
            % dydt flops for computing the slope
        
        % Project k1 to get an estimate for the value of the function
        % at the midpoint of the step,
        ymid = y(i,:) + k1.*hh./2;        
        % and calculate the slope at this midpoint.
        k2 = dydt(tt+hh/2,ymid,varargin{:});
            % 3m flops for projecting
            % 3 flops for setting the time value
            % dydt flops for computing the slope
        
        % Using the new slope k2, 
        % recalculate the projection at the midpoint
        ymid = y(i,:) + k2.*hh./2;        
        % and find another new slope k3.
        k3 = dydt(tt+hh/2,ymid,varargin{:});
            % 3m flops for projecting
            % 3 flops for setting the time value
            % dydt flops for computing the slope
        
        % Using the last slope k3,
        % Re-project all the way to the end of the step
        yend = y(i,:) + k3.*hh;
        % and approximate the slope there.
        k4 = dydt(tt+hh,yend,varargin{:});
            % 3m flops for projecting
            % 2 flops for setting the time value
            % dydt flops for computing the slope
        
        % Use the weighted RK4 formula to refine the approximation.
        phi = (k1 +2*(k2+k3) +k4)/6;
            % 5m flops for the final slope approximation
        
        % With the final slope approximation, project one more time 
        % to find the value of the function for the start of the next step
        y(i+1,:) = y(i,:) +phi*hh;
            % 2m flops for projecting to the next step
            
%--------------------------------------------------------------------------
%% Firing Threshold:
        
        % Check each neuron to see which ones have fired.
        for j=1:m
        
            iNext = i+1;
            
            % The neuron isn't firing at rest
            threshT(iNext,j) = false;
            fireT(iNext,j) = false;
            
            % If we are post fire on a previous neuron,
            if ( tt < resetTime(j) ) && thresholdMask(j)
                % set firing voltage and move to the next neuron.
                y(iNext,j) = holdValue;
                threshT(iNext,j) = true;
                continue
            end
        
            % Any non-post fire neurons we check for threshold potential.
            if ( y(i+1,j) >= threshold )  && thresholdMask(j)
                y(iNext,j) = holdValue;
                threshT(iNext,j) = true;
                fireT(iNext,j) = true;
                resetTime(j) = tt + holdTime;
                continue
            end
        end
%--------------------------------------------------------------------------

        % Set up the next intermediate step in the approximation,
        % and break once we reach the next given time interval
        tt = tt+hh;
        i = i+1;
        if tt>=tend,break,end
            % 1 flop for moving to the next step
            
    end
    
    % Once we have reached the end of each time interval,
    % we record the values for output, and move to the next one.
    np = np+1;
    tOut(np) = tt;
    yOut(np,:) = y(i,:);

    % Once we hit the end of tspan, we break the loop.
    if tt>=tf,break,end
end