%% LIF_net_constBell:
%   Uses leaky integrate and fire, row to column,
%   to model a single, disconnected neuron,
%   with no input,
%   write the solution to a csv file as
%   [ time (s), transmembrane potential (mV) ],
%   with each entry in a separate row and ordered by time,
%   and then display the results.
%##########################################################################
%% Parameters:
%   N = number of neurons in population
%   EI_ratio = ratio for number of excitatory to inhibitory neurons
%   IFS_ratio = ratio for the number of fast to slow inhibitory neurons
%   ----
%   t_init = initial time of desired solution in seconds (s)
%   t_fin = final time of desired solution in seconds (s)
%   step = step size between time points in seconds (s)
%   ----
%   V_init = [ V_init_1, V_init_2, ..., V_init_n ]
%   V_init _i = transmembrane potential of neuron i
%               in milliVolts (mV) at time t_init
%   S_init = [ V_init_1, V_init_2, ..., V_init_n ]
%   S_init _i = transmembrane potential of neuron i
%               in milliVolts (mV) at time t_init 
%   ----
%   In(t) = [ In_1(t); In_2(t); ...; In_n(t) ]
%   In(t) = input control current across membrane 
%           in milliAmps per square centimeter (mA/cm2)
%##########################################################################
%% Pseudocode:
%   ####
%   Parameters:
%   ====
%   Setup:
%   ----
%   Solve:
%   ----
%   Display:
%   ----
%   Save:
%   ####
%##########################################################################
%% Parameters:
N = 3;
EI_ratio = 0;
IFS_ratio = nan;
%--------------------------------------------------------------------------
t_init = 0;     % initial time (s)
t_fin = 5;      % final time (s)
step = 1e-1;    % time interval step size in seconds (s)
%--------------------------------------------------------------------------
in_mean = -60; % transmembrane potential (mV) at time t_init
in_stddev = 0;
%==========================================================================
%% Setup:

% Make parameters for the network of neurons
[ gl, G, El, E, C, ...
    delay, rise, decay, ...
    thresholdMask, threshBases, threshJumps, threshDecays, ...
    refractoryValues, refractoryTimes ] = lif_I1_net();

% Set all constants
dVS =@(t,VS,lastFire) ...
    lif_net( t, VS, lastFire, In, ...
            gl, G, El, E, C, ...
            delay, rise, decay );

%--------------------------------------------------------------------------
% Initial volatges and gate values
V_Init = El;
S_Init = zeros(1,N);  

% Interleave initial variables
VS_Init = zeros(1,2*m);
mask = zeros(1,2*m);
for i=1:m
    VS_Init((2*i)-1) = V_Init(i);
    VS_Init(2*i) = S_Init(i);
end

%==========================================================================
%% Solve:

[ tOut, VOut, fireMask, threshMask ] = ...
    eul_resThresh( dVS, [ t_init, t_fin ], VS_Init, step, ...
                thresholdMask, threshBases, threshJumps, threshDecays, ...
                refractoryValues, refractoryTimes );

%--------------------------------------------------------------------------
%% Display:

% Set up a grid of firing times to pass to the plotter.
firingTimes = nan(length(tOut),2*m);
for j=1:2*m
    subPass = tOut;
    subPass(~fireMask(:,j)) = nan;
    firingTimes(:,j) = subPass';
end

%display(firingTimes);
vid = phasePlot(tOut,VOut,dVS,firingTimes);

% Draw the corresponding firing plot
figure
grid on
hold on
for j=1:m
    fire = j*fireMask(:,(2*j)-1);
    ind = find(fire);
    tFire = tOut(ind);
    iFire = fire(ind);
    firePlot = plot(tFire,iFire,'ro');
end
xlim manual;
xlim([t_init t_fin]);

%--------------------------------------------------------------------------
%% Save:

% Create the names for the save files
firingPlotFileName = sprintf(...
    'firingPlot-LIF_net-noInput_m%d.jpg',m);
phasePlotFileName = sprintf(...
    'phasePlot-LIF_net-noInput_m%d.avi',m);

% Write the firing plot to an image
saveas(gcf,firingPlotFileName);

% and the phase plot to a video.
v = VideoWriter(phasePlotFileName);
open(v);
writeVideo(v,vid);
close(v);

%##########################################################################