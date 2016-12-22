function [ tp, yp, frames ] = eulExODEsysVid( dydt, tspan, y0, h, varargin )
%% eulExODEsysD: solves a system first order ODEs using Euler's method
%   [ tp, yp, frames ] = eulExODEsysD( dydt,tspan,y0,h,varargin )
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
tp = (ti:h:tf)'; 
n = length(tp);
m = length(y0); % Numer of independent variables


% if necessary, add an additional value of t
% so the range goes from t = ti to tf
if tp(n)<tf
    tp(n+1) = tf;
    n = n+1;
end

%==========================================================================
%% Main Algorithm:

yp = ones(n,m);
yp(1,:) = y0';    % preallocate y to improve efficiency
for i = 1:(n-1)         % Step through Euler's method
    yp(i+1,:) = yp(i,:) + ((tp(i+1)-tp(i))*dydt(tp(i),yp(i,:),varargin{:}))';
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