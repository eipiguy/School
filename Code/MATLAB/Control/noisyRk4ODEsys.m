function [ tp, xp ] = noisyRk4ODEsysDisp( dy, tspan, y0, h, e, varargin )
%% noise: noise simulator using rk4 approximation
%   [] = noise( dx, x0, e )
%       Uses the 4th order Runga-Kutta method to approximate
%       the path of the given point withing a noisy system
% input:
%   dx = ODE with independent(t) & dependent(x) inputs
%   tspan = [ti,tf]
%       ti = initial time
%       tf = final time
%   OR tspan = [ti,t1,t2...,tf] points to approximate the function at
%   x0 = initial values of dependent variables
%   h = step size
%   p1,p2,... = additional parameters used by dx
%   e = error/noise function to be added into differential equation
% output:
%   tp = vector of independent variables
%   xp = vector of solution for dependent variables
%% Psuedo-Code:
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
if nargin<5
    error('ODE, time span, initial conditions, and step size, and noise function required.');
end
% Make sure initial and final times are in increasing order
if any(diff(tspan)<=0),error('tspan not in ascending order');end

%==========================================================================
%% Variable Declarations:

% Add in noise to dynamic system
    function [ dxdt ] = dx( ty, yi, varargin )
        dxdt = dy(ty,yi,varargin{:}) + e(ty,yi,varargin{:}) ;
    end

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

% Set up the other initial parameters for the method
tt = ti; x(1,:) = y0;
np = 1; tp(np) = tt; xp(np,:) = x(1,:);
i = 1;

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
        k1 = dx(tt,x(i,:),varargin{:})';
            % dydt flops for computing the slope
        
        % Project k1 to get an estimate for the value of the function
        % at the midpoint of the step,
        xmid = x(i,:) + k1.*hh./2;        
        % and calculate the slope at this midpoint.
        k2 = dx(tt+hh/2,xmid,varargin{:})';
            % 3m flops for projecting
            % 3 flops for setting the time value
            % dydt flops for computing the slope
        
        % Using the new slope k2, 
        % recalculate the projection at the midpoint
        xmid = x(i,:) + k2.*hh./2;        
        % and find another new slope k3.
        k3 = dx(tt+hh/2,xmid,varargin{:})';
            % 3m flops for projecting
            % 3 flops for setting the time value
            % dydt flops for computing the slope
        
        % Using the last slope k3,
        % Re-project all the way to the end of the step
        xend = x(i,:) + k3.*hh;
        % and approximate the slope there.
        k4 = dx(tt+hh,xend,varargin{:})';
            % 3m flops for projecting
            % 2 flops for setting the time value
            % dydt flops for computing the slope
        
        % Use the weighted RK4 formula to refine the approximation.
        phi = (k1 +2*(k2+k3) +k4)/6;
            % 5m flops for the final slope approximation
        
        % With the final slope approximation, project one more time 
        % to find the value of the function for the start of the next step
        x(i+1,:) = x(i,:) +phi*hh;
            % 2m flops for projecting to the next step
        
        % Move to the next intermediate step in the approximation,
        % and break once we reach the next given time interval
        tt = tt+hh;
        i = i+1;
        if tt>=tend,break,end
            % 1 flop for moving to the next step
    end
    
    % Once we have reached the end of each time interval,
    % we record the values for output, and move to the next one.
    np = np+1; tp(np) = tt; xp(np,:) = x(i,:);
    
    % Once we hit the end of tspan, we break the loop.
    if tt>=tf,break,end
end
%##########################################################################
end