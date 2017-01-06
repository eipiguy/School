function [ tOut,yOut,threshOut,fireMask ] = ...
    eul_refThresh( dydt,tspan,y0,h,thresholdMask,...
                    threshBases,threshJumps,threshDecays, ...
                    refractoryValues,refractoryTimes,varargin )
%% eul_refThresh: solves a system of ODEs with a thresheld refractory
%   [ tOut,yOut,threshOut,fireMask ] = eul_refThresh( dydt,tspan,y0,h, ...
%                   thresholdMask,threshBase,threshJump,threshDecay,
%                   refractoryPeriod,holdTime,varargin )
%       Uses the implicit euler method
%       to numerically solve a system of first order ODEs,
%       a number of which change in response to a threshold,
%       as indicated by the thresholdMask. Passing a threshold
%       produces an increase in the threshold of that variable,
%       which decays back to its normal value exponentially.
%       In addition, passing a threshold sets that variable
%       to a given refractory value for a period of time,
%       before giving control back over to the solver.
%==========================================================================
%% Input:
%   dydt = differential equation with independent(t) & dependent(y) inputs
%   tspan = [ti,tf]
%       ti = initial time
%       tf = final time
%   OR tspan = [ti,t1,t2...,tf] points to approximate the function at
%   y0 = initial values of dependent variables
%   h = step size
%--------------------------------------------------------------------------
%   thresholdMask = row vector of zeros and ones 
%                   where a one indicates to threshold check
%                   the corresponding variable in the system
%   threshBases = vector of base threshold values
%   threshJumps = vector of values to increase threshold after fire
%--------------------------------------------------------------------------
%   refractoryValues = vector of firing refractory hold values
%   refractoryTimes = vector of times to hold for each variable
%--------------------------------------------------------------------------
%   p1,p2,... = additional parameters used by dydt
%==========================================================================
%% Output:
%   tOut = vector of independent time values for the matching solution
%   yOut = vector of phase solution for dependent variables
%   threshOut = vector of threshold values for each variable at each time
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
if nargin<10
    error('Enter all input parameters.');
end
% Make sure initial and final times are in increasing order
if any(diff(tspan)<=0),error('tspan not in ascending order');end

m = length(y0);                 % Number of variables in the system

% Make sure the inputs have the correct number of variables
if length(thresholdMask) ~= m
    error('y0 and thresholdMask must have the same number of variables');
end
if length(refractoryValues) ~= m
    error('y0 and refractoryValues must have the same number of variables');
end
if length(refractoryTimes) ~= m
    error('y0 and refractoryTimes must have the same number of variables');
end

%==========================================================================
%% Variable Declarations:

% Setting up the steps for the time span:
n = length(tspan);              % Number of steps between endpoints
ti = tspan(1); tf = tspan(n);   % Set variables first time and last time

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

%--------------------------------------------------------------------------
% Set up the other initial parameters for the method:

% Inner solver starting conditions
i=1;            % curent partition number
tt = ti;        % current time
y(1,:) = y0;    % current state vector

% Threshold condition initializations
resetTime = zeros(m);           % time before refractory over
fireMask = zeros(n,m);

lastFire = nan(1,m);
curThresh = threshBases;
fireStack = zeros(1,2*m);

%--------------------------------------------------------------------------
% Initial Threshold Check
for j=1:m
    % If we start above a threshold value, 
    % then we change the initial value to match firing behavior    
    if (y0(j)>=curThresh(j)) && thresholdMask(j)==1;
                
        % Set the firing flags, and update the reset time
        fireMask(1,j) = true;
        resetTime(j) = tt + refractoryTimes(j);
        lastFire(j) = tt;
        fireStack(j) = fireStack(j)+1;
            % 1 flop for changing the reset time
        
        % Set the refractory value, and change the threshold
        y(1,j) = refractoryValues(j);
        curThresh(j) = curThresh(j) + threshJumps(j);
        fireStack(j) = threshJumps(j);
            % 1 flop for changing the threshold
    end
end
%--------------------------------------------------------------------------

% Set output to partition starting conditions
np = 1; tOut(np) = tt; yOut(np,:) = y(1,:); threshOut(np,:) = curThresh;

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
            
        iNext = i+1;    % next partition number
        
        % Use Euler's method to project to the end of the time step
        y(iNext,:) = y(i,:) + ( hh*dydt(tt,y(i,:),lastFire,varargin{:}) );
            % dydt flops for computing the slope
            % 2m flops for multiplying by step size and adding
            
%--------------------------------------------------------------------------
%% Firing Threshold:

        % Check each neuron to see which ones have fired.
        for j=1:m
            
            % The threshold decays after every fire
            if ~isnan(lastFire(j)) && thresholdMask(j)==1
                curThresh(j) = threshBases(j)+(fireStack(j)*exp(-(tt-lastFire(j))/threshDecays(j)));
            end
            
            % If we are post fire on a previous neuron,
            if ( tt < resetTime(j) ) && thresholdMask(j)==1
                % set firing voltage and move to the next neuron.
                y(iNext,j) = refractoryValues(j); 
                continue
            end
        
            % Any non-post fire neurons we check for threshold potential.
            if ( y(iNext,j) >= curThresh(j) )  && thresholdMask(j)==1
                y(iNext,j) = refractoryValues(j);
                curThresh(j) = curThresh(j) + threshJumps(j);
                fireStack(j) = curThresh(j) - threshBases(j);
                fireMask(np+1,j) = true;
                resetTime(j) = tt + refractoryTimes(j);
                lastFire(j) = tt;
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
    threshOut(np,:) = curThresh;

    % Once we hit the end of tspan, we break the loop.
    if tt>=tf,break,end
end