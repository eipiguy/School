function redraw(frame)
%% redraw: function that runs for each frame in videofig.m
%   redraw(frame)
%       This is run each time a frame is rendered using videofig.
%       This particular setup is used to animate the thetaNeuron model
%       on top of the unit circle.
% input:
%   frame = current frame number
%% Pseudo-Code:
%   ####
%   Main Algorithm:
%   ####
%##########################################################################
%% Main Algorithm:

    % start with a fresh background circle each frame.
    hold off;
    unitCircle(0.01);
    xlim([ -1.5 1.5 ]);
    ylim([ -1.5 1.5 ]);
    
    % approximate the ODE solution and store it in persistent variables
    persistent y;
    persistent y1;
    persistent y2;
    if length(y1)<2
        [t,y] = rk4ODEsys(@thetaNeuron,[0,25],0,0.05,@currentIn);
        y1 = cos(y);
        y2 = sin(y);
    end
    
    % save the background image and place the current point on the circle
    hold on;
    plot(y1(frame),y2(frame),'ro');
    
%##########################################################################
end