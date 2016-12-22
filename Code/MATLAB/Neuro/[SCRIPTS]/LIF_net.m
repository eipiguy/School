%% LIF_net_zeroIn:
%   Uses leaky integrate and fire, row to column,
%   to model a single, disconnected neuron,
%   with no input,
%   write the solution to a csv file as
%   [ time (s), transmembrane potential (mV) ],
%   with each entry in a separate row and ordered by time,
%   and then display the results.
%##########################################################################
%% Parameters:
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
%   ----
%   g = [ g1; g2; ...; gn ]
%   gi = membrane conductance in miliSiemens per centimeter^2 (mS/cm2)
%   E = [ E1; E2; ...; En ]
%   Ei = leak equilibrium potential in miliVolts (mV)
%   C = [ C1; C2; ...; Cn ]
%   Ci = membrane capacitence in miliFarads per centimeter^2 (mF/cm2)
%   ----
%   thresh = [ thresh1; thresh2; thresh3; ...; threshn ]
%   thresh = threshold of transmembrane potential in miliVolts (mV)
%               that produces firing behavior
%   holdV = [ holdV1; holdV2; ...; holdVn ]
%   holdV = transmembrane potential in miliVolts (mV) 
%               after refractory period
%   holdT = [ holdT1; holdT2; ...; holdTn ]
%   holdT = time of firing refractory period in seconds(s)
%   ----
%   delay = [ delay1; delay2; ...; delayn ];
%   rise = [ rise1; rise2; ...; risen ];
%   decay = [ deay2; decay2; ...; decayn ];
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
t_init = 0;     % initial time (s)
t_fin = 5;      % final time (s)
step = 1e-1;    % time interval step size in seconds (s)
%--------------------------------------------------------------------------
V_Init = -60*ones(1,3);  % transmembrane potential (mV) at time t_init
S_Init = zeros(1,3);      % Initial gate variable values at time t_init
%--------------------------------------------------------------------------
cIn = 3.5*ones(3,1);           % constant input current vector
cIn(2:3) = 0;
In =@(t) constIn(t,cIn);    % input control current (mA/cm2)
%--------------------------------------------------------------------------
gl = 0.25*ones(3,1);    % membrane conductivity (mS/cm2)
G = [   0,  0,  0;
        0.1,  0,  0;
        0,  0.1,  0   ];  % maximum net-synapse conductance matrix (mS/cm2)
El = -68.6*ones(3,1);   % transmembrane potential (mV) at equilibrium
E = 0*ones(3,1);           % synaptic equilibrium voltage (mV)
C = 0.1*ones(3,1);      % membrane capacitance (mF/cm2)
%--------------------------------------------------------------------------
thresh = -55*ones(6,1);  % threshold potential (mV)
holdV = -68.6*ones(6,1);    % firing refractory hold potential (mV)
holdT = 1.5e-3*ones(6,1);   % time of post fire refractory period (s)
%--------------------------------------------------------------------------
delay = 0.2;
rise = 0.2;
decay = 0.75;
%==========================================================================
%% Setup:

m = length(V_Init);

% Set all constants
dVS =@(t,VS,lastFire) lif_net( t,VS,lastFire,In,delay,rise,decay,gl,G,El,E,C );

% Interleave initial variables
VS_Init = zeros(1,2*m);
mask = zeros(1,2*m);
for i=1:m
    VS_Init((2*i)-1) = V_Init(i);
    VS_Init(2*i) = S_Init(i);
    % Thresholds should only be checked on voltage variables,
    % not on the gate variables
    mask((2*i)-1) = true;
end


%--------------------------------------------------------------------------
%% Solve:

[ tOut, VOut, fireMask, threshMask ] = eul_resThresh( ...
                                dVS, [ t_init, t_fin ], VS_Init, ...
                                step, mask, thresh, holdV, holdT );

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