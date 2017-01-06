function [ frames ] = phaseInPlot( tSol, In, ySol, threshSol, dydt, varargin )
%% Template: One line summary of the function goes here
%   Detailed explanation goes here
%##########################################################################
%% Input:
%   tSol = column vector of independent time variable in given units
%   ySol = row vectors of dependent variables 
%           in a column corresponding to each time in tSol
%==========================================================================
%% OutSolut:
%   y = describe each outSolut variable in units (unit abbreviation)
%##########################################################################
%% Pseudocode:
%   ####
%   Input Check:
%   ====
%   Variable Declaration:
%   ----
%   Figure Setup:
%   ====
%   Time Plots:
%   ----
%   Phase Planes:
%   ====
%   Video:
%##########################################################################
%% Input Check:

%==========================================================================
%% Variable Declaration:

tPTS = 10;  % Number of quivers on the time axis
yPTS = 10;  % Number of quivers on the solution axis

m = size(ySol,2);           % Number of y components in solution
ySpan = zeros(2,m);         % Initialize min/max value storage

ti = tSol(1);               % Initial solution time
tf = tSol(length(tSol));    % Final solution time
tL = length(tSol);          % Number of time samples in solution

tPPDiv = floor(tL/tPTS);    % Number of time points between grid points
tDivInd = 1:tPPDiv:tL;      % Solution time index for each grid point
tRV = tSol(tDivInd);        % Sample time point values in a row vector
tRL = length(tRV);          % Number of resulting "key" times

%--------------------------------------------------------------------------
%% Figure Setup:

% Clear all current figures and turn off hold
cla;
clf;
hold on;

% Get the current screen size
figHist = gcf;
scrsz = get(groot,'ScreenSize');

% Resize the figure window in the middle of the screen at half it's size
set(figHist,'Position', ...
    [0 0 scrsz(3) scrsz(4)]);

% Set up the grid of sub plots for inputs, time, and phase.
for i=1:(m^(2)+m)
    sub(i) = subplot(m,m+1,i);
    set(sub(i),'Visible','off');
end

%==========================================================================
%% Input Plots:
    
% Place each input on the leftmost column
for i=1:m
    if mod(i,2)==1
        subplot(sub(1+(i-1)*(m+1)));    % Point to the subplot
        % Draw the variable solution and threshold vs time
        plot(tSol,In(((i+1)./2),:));
        ylim([ min([0,min(In(((i+1)./2),:))]), max([0,max(In(((i+1)./2),:))]) ]);
        xLab = xlabel('t');                 % Label the x-axis as time
        if mod(i,2)==1
            yLab = ylabel(sprintf('In(%d)',(i+1)./2));  % Label the y-axis with the index
        end
        set(gca,'FontSize',4+(8/m));        % Size scales with number of plots
        xLab.FontSize = 4+(8/m);
        yLab.FontSize = 4+(8/m);
    end
    
    % Turn on hold so the solutions stay visible
    %xlim manual;
    %ylim manual;
    hold(sub(1+(i-1)*(m+1)),'on');
    grid(sub(1+(i-1)*(m+1)),'on');
end

%==========================================================================
%% Time Plots:

% Place each time plot solution in the second column.
for i=1:m
    subplot(sub(2+(i-1)*(m+1)));    % Point to the subplot
    % Draw the variable solution and threshold vs time
    plot(tSol,ySol(:,i),tSol,threshSol(:,i),'--');
    %plot(tSol,threshSol(:,i),'--');
    xLab = xlabel('t');                 % Label the x-axis as time
    if mod(i,2)==1
        yLab = ylabel(sprintf('V(%d)',(i+1)./2));  % Label the y-axis with the index
    else
        yLab = ylabel(sprintf('S(%d)',i./2));  % Label the y-axis with the index
    end
    set(gca,'FontSize',4+(8/m));        % Size scales with number of plots
    xLab.FontSize = 4+(8/m);
    yLab.FontSize = 4+(8/m);
    
    
    % Set the min/max values for each solution curve
    ySpan(:,i) = ylim';
    
    % Turn on hold so the solutions stay visible
    xlim manual;
    ylim manual;
    hold(sub(2+(i-1)*(m+1)),'on');
    grid(sub(2+(i-1)*(m+1)),'on');
end

% Set the difference between solution grid points
% after the min and max have been found
for i=1:m
    ySpan(3,i) = (ySpan(2,i)-ySpan(1,i))/yPTS;
end

%--------------------------------------------------------------------------
%% Phase Planes:

if m >= 2
% Place phase graphs next to their corresponding time solutions
for i=1:m-1
    for j=i+1:m
        subplot(sub(((i-1)*(m+1))+j+1));  % Point to the subplot
        plot(ySol(:,j),ySol(:,i));  % Plot the orbit in the phase plane
        if mod(j,2)==1
            xLab = xlabel(sprintf('V(%d)',(j+1)./2));  % Label the y-axis with the index
        else
            xLab = xlabel(sprintf('S(%d)',j./2));  % Label the y-axis with the index
        end
        if mod(i,2)==1
            yLab = ylabel(sprintf('V(%d)',(i+1)./2));  % Label the y-axis with the index
        else
            yLab = ylabel(sprintf('S(%d)',i./2));  % Label the y-axis with the index
        end
        set(gca,'FontSize',4+(8/m));
        xLab.FontSize = 4+(8/m);
        yLab.FontSize = 4+(8/m);
        
        % Turn on hold so the solutions stay visible
        xlim manual;
        ylim manual;
        hold(sub(((i-1)*(m+1))+j+1),'on');
        grid(sub(((i-1)*(m+1))+j+1),'on');
    end
end
end

%==========================================================================
%% Make Video Frames:

% For each component of the solution, 
% make a video of the background quiver plot
for i=1:m
    
    counter = 0;
    
    % Point at the correct subplot, 
    subplot(sub(2+(i-1)*(m+1)));
    
    % For each component, start by finding the static quiver plot
    % of that component's solution with respect to time
    for t=tDivInd
        counter = counter+1;
            
        % Set the row vector of "key" grid points 
        % for the solution values
        yRV = ySpan(1,i):ySpan(3,i):ySpan(2,i);
            
        % If the quiver data variable doesn't exist, initialize it
        if exist('dGrid','var')==0
            dGrid = zeros(length(tRV),length(yRV));
        end
            
        % For each grid point in the row vector,
        for yIndex=1:length(yRV)
                
            % Create the y input to feed into dydt                
            yIn = ySol(t,:);        % Start with the current solution
            yIn(i) = yRV(yIndex);   % Set the current component to 
                                    % the solution at that grid point
                
            
            % Compute the associated slope and set the place in the grid
            dY = dydt( tSol(t), yIn, varargin{:} );
            dGrid(counter,yIndex,i) = dY(i);
        end
    end
    
    dispGrid = dGrid(:,:,i)';
    
    % Set up the vertices of the quivers
    [ inMeshX, inMeshY ] = meshgrid(tRV,yRV);
    
    q = quiver(inMeshX,inMeshY,ones(size(dispGrid)),dispGrid);
    q.Color = 'k';
    q.ShowArrowHead = 'off';
    q.AlignVertexCenters = 'on';
end

% For each time in the solution,
% make a 
for t=1:tL;
    % Place each time plot solution in the second column.        
    for i=1:m  
        % if we are at an odd number, plot the input as well
        if mod(i,2)==1
            subplot(sub(1+(i-1)*(m+1)));
            plot(tSol(t),In((i+1)./2,t),'ro')
        end
        
        % Point at the correct subplot, 
        subplot(sub(2+(i-1)*(m+1)));
            
        % Overlay the current time point
        plot(tSol(t),ySol(t,i),'ro');
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
                subplot(sub(((i-1)*(m+1))+j+1));
            
                % Create the row vectors of key grid points for each axis
                xPhRV = ySpan(1,j):ySpan(3,j):ySpan(2,j);
                ySolhRV = ySpan(1,i):ySpan(3,i):ySpan(2,i);
            
                % Record the number of key grid points on each axis
                xRVL = length(xPhRV);
                yRVL = length(ySolhRV);
            
                % Compute the current quiver grid
                dPhGridX = zeros(xRVL,yRVL);
                dPhGridY = zeros(xRVL,yRVL);
                for xIndex=1:xRVL
                    for yIndex=1:yRVL
                
                        % Create the y input to feed into dydt
                        % With each input slot as a row
                
                        yIn = ySol(t,:);    % Start with current y values
                        
                        % Set the current component to match a grid point
                        yIn(j) = xPhRV(xIndex);
                        yIn(i) = ySolhRV(yIndex);
                
                        % Compute the associated slope 
                        % and set the place in the grid
                        dY = dydt( tSol(t), yIn, varargin{:} );
                        dPhGridX(xIndex,yIndex) = dY(j);
                        dPhGridY(xIndex,yIndex) = dY(i);
                    end
                end
                
                dispGridX = dPhGridX';
                dispGridY = dPhGridY';
            
                % Set up the vertices of the quivers
                [ inMeshX, inMeshY ] = meshgrid(xPhRV,ySolhRV);
                
                % Draw the current quiver grid
                dq(i,j) = quiver(inMeshX,inMeshY,dispGridX,dispGridY);
                xlim([ySpan(1,j),ySpan(2,j)]);
                ylim([ySpan(1,i),ySpan(2,i)]);
                dq(i,j).Color = 'k';
                dq(i,j).ShowArrowHead = 'off';
                dq(i,j).AlignVertexCenters = 'on';
        
                % Overlay the current time point
                plot(ySol(t,j),ySol(t,i),'ro');
            end
        end
    end
    % Record the current frame
    frames(t) = getframe(figHist);
end

%##########################################################################
end

