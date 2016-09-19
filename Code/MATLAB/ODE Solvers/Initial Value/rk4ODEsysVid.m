function [ tp,yp,frames ] = rk4ODEsysVid( dydt,tspan,y0,h,varargin )
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
        
        % Move to the next intermediate step in the approximation,
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
    if (floor(i./m) > mod(i,m)+1)
        set(sub(i),'Visible','off');
    end
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
            
            ySpan(3,i) = (ySpan(2,i)-ySpan(1,i))/10;
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

    if m >= 2
    % Place phase graphs next to their corresponding time solutions
    for i=1:m-1
        for j=i+1:m
            % Start with the background dynamics
            
            % Point at the correct subplot, 
            subplot(sub(((i-1)*m)+j));
            
            xPhDiv = (ySpan(2,j)-ySpan(1,j))/10;
            xPhRV = ySpan(1,j):xPhDiv:ySpan(2,j);
            yPhDiv = (ySpan(2,i)-ySpan(1,i))/10;
            yPhRV = ySpan(1,i):yPhDiv:ySpan(2,i);
            
            % Compute the matrix of the current dY with
            % time slots stored as rows again
            dPhGridX = zeros(length(xPhRV),length(yPhRV));
            dPhGridY = zeros(length(xPhRV),length(yPhRV));
            for xM=1:length(xPhRV)
                for yM=1:length(yPhRV)
                
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
            
            dq.Visible = 'off';
            dq = quiver(inMeshX,inMeshY,dispGridX,dispGridY);
            dq.Color = 'k';
            dq.Marker = '+';
            dq.MarkerSize = 3;
            dq.ShowArrowHead = 'off';
            %q.MaxHeadSize = min(tDiv,ySpan(3,i));
            %dq.AutoScale = 'off';
            dq.AlignVertexCenters = 'on';
        
            % Overlay the current time point
            plot(yp(t,j),yp(t,i),'ro');
        end
    end
    end
    frames(t) = getframe(figHist);
end
%##########################################################################
end