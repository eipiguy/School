function [ tp,yp,fireT,threshT,frames ] = rk4NeuroSysVid( dydt,tspan,y0,h,threshold,fireV,holdTime,varargin )
%% rk4ODEsysVid: solves a system of first order ODEs, and plots the result
%   [ t,y ] = rk4sys( dydt,tspan,y0,h,varargin )
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
%   threshold = threshold potential for each neuron
%   fireV = voltage to hold each neuron to simulate firing
%   p1,p2,... = additional parameters used by dydt
% output:
%   tp = vector of independent variables
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

% Set up the other initial parameters for the method
tt = ti; y(1,:) = y0;
np = 1; tp(np) = tt; yp(np,:) = y(1,:);
i = 1; 
resetTime = zeros(m);
threshCurrent = threshold;
threshT = zeros(n,m);
fireT = zeros(n,m);
threshIT = zeros(n,m);
fireIT = zeros(n,m);

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
        k1 = dydt(tt,y(i,:),varargin{:})';
            % dydt flops for computing the slope
        
        % Project k1 to get an estimate for the value of the function
        % at the midpoint of the step,
        ymid = y(i,:) + k1.*hh./2;        
        % and calculate the slope at this midpoint.
        k2 = dydt(tt+hh/2,ymid,varargin{:})';
            % 3m flops for projecting
            % 3 flops for setting the time value
            % dydt flops for computing the slope
        
        % Using the new slope k2, 
        % recalculate the projection at the midpoint
        ymid = y(i,:) + k2.*hh./2;        
        % and find another new slope k3.
        k3 = dydt(tt+hh/2,ymid,varargin{:})';
            % 3m flops for projecting
            % 3 flops for setting the time value
            % dydt flops for computing the slope
        
        % Using the last slope k3,
        % Re-project all the way to the end of the step
        yend = y(i,:) + k3.*hh;
        % and approximate the slope there.
        k4 = dydt(tt+hh,yend,varargin{:})';
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
%%Intermediate Threshold
        
        % Check each neuron to see which ones have fired.
        for j=1:m
        
            % The neuron isn't firing at rest
            threshIT(i,j) = false;
            fireIT(i,j) = false;
            
            % If we are post fire on a previous neuron,
            if ( tt < resetTime(j) )
                % set firing voltage and move to the next neuron.
                y(i+1,j) = fireV;
                threshIT(i,j) = true;
                continue
            end
        
            % Any non-post fire neurons we check for threshold potential.
            if ( y(i+1,j) >= threshCurrent )
                y(i+1,j) = fireV;
                threshIT(i,j) = true;
                fireIT(i,j) = true;
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
    np = np+1; tp(np) = tt; yp(np,:) = y(i,:);

    % Once we hit the end of tspan, we break the loop.
    if tt>=tf,break,end

%--------------------------------------------------------------------------
%%Firing Threshold
        
        % Check each neuron to see which ones have fired.
        for j=1:m
        
            % The neuron isn't firing at rest
            threshT(np,j) = false;
            fireT(np,j) = false;
            
            % If we are post fire on a previous neuron,
            if ( tt < resetTime(j) )
                % keep firing voltage and move to the next neuron.
                yp(np,j) = fireV;
                threshT(np,j) = true;
                fireT(np,j) = true;
                continue
            end
        
            % Any non-post fire neurons we check for threshold potential.
            if ( yp(np,j) >= threshCurrent )
                yp(np,j) = fireV;
                resetTime(j) = tt + holdTime;
                threshT(np,j) = true;
                continue
            end
        end
%--------------------------------------------------------------------------
end

%==========================================================================
%% Plot Results:

% Clear all current figures and turn off hold
cla;
clf;
hold on;

% get the current screen size, and resize the figure window
% to fit directly in the middle of the screen at half it's size

figHist = gcf;
scrsz = get(groot,'ScreenSize');
set(figHist,'Position',[scrsz(3)/6 scrsz(4)/6 (2*scrsz(3))/3 (2*scrsz(4))/3]);

% Set up the grid of sub plots for each of the time and phase plots.
for i=1:(m^(2))
    sub(i) = subplot(m,m,i);
    %if (floor(i./m) > mod(i,m)+1)
        set(sub(i),'Visible','off');
    %end
end

ySpan = zeros(2,m);
% Place each time plot solution in the left column.
for i=1:m
    subplot(sub(1+(i-1)*m));
    plot(tp,yp(:,i));
    xlabel('t');
    ylabel(sprintf('y(%d)',i));
    
    % Set the span of the y-axis to be the 
    % column with the corresponding index
    ySpan(:,i) = ylim';
    
    %display(ySpan);
    xlim manual;
    ylim manual;
    hold(sub(1+(i-1)*m),'on');
end
for i=1:m
    ySpan(3,i) = (ySpan(2,i)-ySpan(1,i))/10;
end
%display(ySpan);

if m >= 2
% Place phase graphs next to their corresponding time solutions
for i=1:m-1
    for j=i+1:m
        subplot(sub(((i-1)*m)+j));
        plot(yp(:,j),yp(:,i));
        xlabel(sprintf('y(%d)',j));
        ylabel(sprintf('y(%d)',i));
        xlim manual;
        ylim manual;
        hold(sub(((i-1)*m)+j),'on');        
    end
end
end
%==========================================================================
%% Make Video Frames:

% Set the quivers' t and y seperation
tDiv = (tf-ti)/20;

% Create the row vectors for the sample points
% of the time and the current y component
tRV = ti:tDiv:tf;
%display(tRV);
cDiv = floor(length(tp)/10);
counter = 0;

for t=1:length(tp);
    % Place each time plot solution in the left column.        
    for i=1:m
        
        % Point at the correct subplot, 
        subplot(sub(1+(i-1)*m));
        
        % If we are at a "key" time,
        % then compute a new column of dY
        if (t==1)||(mod(t,cDiv)==0)
            counter = counter+1;
            
            yRV = ySpan(1,i):ySpan(3,i):ySpan(2,i);
            
            % Compute the matrix of the current dY with
            % time slots stored as rows again
            if exist('dGrid','var')==0
                dGrid = zeros(length(tRV),length(yRV));
            end
            for yM=1:length(yRV)
                
                % Create the y input to feed into dydt
                % With each input slot as a row
                
                yIn = yp(t,:);      % Start with the current y values
                yIn(i) = yRV(yM);   % Set the current component to 
                                    % match a grid point
                %display(yIn);
                
                % Compute the associated slope 
                % and set the place in the grid
                dY = dydt( tp(t), yIn, varargin{:} );
                dGrid(counter,yM) = dY(i);
            end
            %display(dGrid);
        
            % Normalize the solution vectior's components
            % to fit in between each quiver
            %dMax = max(max(abs(dGrid)));
            %display(dMax);
            dispGrid = dGrid(counter,:)';%/dMax;
        
            % Set up the vertices of the quivers
            [ inMeshX, inMeshY ] = meshgrid(tp(t),yRV);
            %display(inMeshX);
            %display(inMeshY);
            
            q = quiver(inMeshX,inMeshY,ones(size(dispGrid)),dispGrid);
            q.Color = 'k';
            q.Marker = '+';
            q.MarkerSize = 3;
            q.ShowArrowHead = 'off';
            %q.MaxHeadSize = max(tDiv,ySpan(3,i));
            %q.AutoScale = 'off';
            q.AlignVertexCenters = 'on';
        end
            
        % Overlay the current time point
        plot(tp(t),yp(t,i),'ro');
    end

    % Place phase graphs next to their corresponding time solutions
    if m >= 2
        if exist('dq','var')
            delete(dq);
        end
        for i=1:m-1
            for j=i+1:m
                % Start with the background dynamics
            
                % Point at the correct subplot, 
                subplot(sub(((i-1)*m)+j));
            
                %display(ySpan);
                xPhRV = ySpan(1,j):ySpan(3,j):ySpan(2,j);
                yPhRV = ySpan(1,i):ySpan(3,i):ySpan(2,i);
            
                xRVL = length(xPhRV);
                yRVL = length(yPhRV);
            
                % Compute the matrix of the current dY with
                % time slots stored as rows again
                dPhGridX = zeros(xRVL,yRVL);
                dPhGridY = zeros(xRVL,yRVL);
                for xM=1:xRVL
                    for yM=1:yRVL
                
                        % Create the y input to feed into dydt
                        % With each input slot as a row
                
                        yIn = yp(t,:);      % Start with the current y values
                        % Set the current component to match a grid point
                        yIn(j) = xPhRV(xM);
                        yIn(i) = yPhRV(yM);
                        %display(yIn);
                
                        % Compute the associated slope 
                        % and set the place in the grid
                        dY = dydt( tp(t), yIn, varargin{:} );
                        dPhGridX(xM,yM) = dY(j);
                        dPhGridY(xM,yM) = dY(i);
                    end
                end
                %display(dPhGridX);
                %display(dPhGridY);
        
                % Normalize the solution vectior's components
                % to fit in between each quiver
                %dMaxX = max(max(abs(dPhGridX)));
                %dMaxY = max(max(abs(dPhGridY)));
                %display(dMaxX);
                %display(dMaxY);
                dispGridX = dPhGridX';%/dMaxX;
                dispGridY = dPhGridY';%/dMaxY;
            
                % Set up the vertices of the quivers
                [ inMeshX, inMeshY ] = meshgrid(xPhRV,yPhRV);
                %display(inMeshX);
                %display(inMeshY);
                
                dq(i,j) = quiver(inMeshX,inMeshY,dispGridX,dispGridY);
                dq(i,j).Color = 'k';
                dq(i,j).Marker = '+';
                dq(i,j).MarkerSize = 3;
                dq(i,j).ShowArrowHead = 'off';
                %dq(i,j).MaxHeadSize = min(tDiv,ySpan(3,i));
                %dq(i,j).AutoScale = 'off';
                dq(i,j).AlignVertexCenters = 'on';
        
                % Overlay the current time point
                plot(yp(t,j),yp(t,i),'ro');
            end
        end
    end
    frames(t) = getframe(figHist);
end
%##########################################################################
end