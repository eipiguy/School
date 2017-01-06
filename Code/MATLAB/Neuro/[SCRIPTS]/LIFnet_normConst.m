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
E_percent = 0.4;
%EIcon_percent_matrix = zeros(N);
EIcon_percent_matrix = [    1, 1;
                            1, 1  ];
%--------------------------------------------------------------------------
t_init = 0;     % initial time (ms)
t_fin = 10;      % final time (ms)
step = 0.01;    % time interval step size in miliseconds (ms)
%--------------------------------------------------------------------------
in_mean = 1e-6; % transmembrane potential (mV) at time t_init
in_stddev = 0;
%==========================================================================
%% Setup:

% Make parameters for the network of neurons
[ gl, G, El, E, C, ...
    delay, rise, decay, ...
    thresholdMask, threshBases, threshJumps, threshDecayTConsts, ...
    refractoryValues, refractoryTimes ] = lif_EI1_net(N,E_percent,EIcon_percent_matrix);

% Set up the input
Ins = zeros(N,1);
for i=1:N
    Ins(i) = normrnd(in_mean,in_stddev);
end
In =@(t) constIn(t,Ins); 

% Set all constants
dVS =@(t,VS,fireTimes) ...
    lif_net( t, VS, fireTimes, In, ...
            gl, G, El, E, C, ...
            delay, rise, decay );

%--------------------------------------------------------------------------
% Initial volatges and gate values
V_Init = El;
S_Init = zeros(1,N);  

% Interleave initial variables
VS_Init = zeros(1,2*N);
for i=1:N
    VS_Init((2*i)-1) = V_Init(i);
    VS_Init(2*i) = S_Init(i);
end

%==========================================================================
%% Solve:

[ tOut, VOut, threshOut, fireMask ] = ...
    eul_refThresh( dVS, [ t_init, t_fin ], VS_Init, step, ...
                thresholdMask, threshBases, threshJumps, threshDecayTConsts, ...
                refractoryValues, refractoryTimes );
            
%--------------------------------------------------------------------------
%% Display:

% Set up a grid of firing times to pass to the plotter.
firingTimes = nan(length(tOut),2*N);
for j=1:2*N
    subPass = tOut;
    subPass(~fireMask(:,j)) = nan;
    firingTimes(:,j) = subPass';
end

% Remove nan solutions from the threshold solutions
% (variables without a threshold)
for i=1:2*N
    if ~thresholdMask(i)
        threshOut(:,i) = nan(1,length(tOut));
    end
end

vid = phaseInPlot(tOut,In(tOut),VOut,threshOut,dVS,firingTimes);

% Draw the corresponding firing plot
figure
grid on
hold on
for j=1:N
    fire = j*fireMask(:,(2*j)-1);
    ind = find(fire);
    tFire = tOut(ind);
    iFire = fire(ind);
    firePlot = plot(tFire,iFire,'ro');
    xLab = xlabel(sprintf('time (ms)'));            % label x-axis as time
    yLab = ylabel(sprintf('neuron #',(i+1)/2));     % label the y-axis 
end                                                 %   with the index
xlim manual;
xlim([t_init t_fin]);
firingFigure = gcf;

figure
dGraph = digraph(G);
plot(dGraph);

%--------------------------------------------------------------------------
%% Save:

oldFolder = pwd;
parentName = oldFolder;

childName = sprintf('Neuro_Results');
childName = fullfile(parentName,childName);
childName = strrep(childName,'.','_');
if ~exist(childName,'dir')
    mkdir(childName);
end
cd(childName);
parentName = childName;

childName = sprintf('E_Percent-%0.2f',E_percent);
childName = fullfile(parentName,childName);
childName = strrep(childName,'.','_');
if ~exist(childName,'dir')
    mkdir(childName);
end
cd(childName);
parentName = childName;

childName = sprintf('LIF_net-%d',N);
childName = fullfile(parentName,childName);
if ~exist(childName,'dir')
    mkdir(childName);
end
cd(childName);
parentName = childName;

childName = sprintf('Const_In');
childName = fullfile(parentName,childName);
if ~exist(childName,'dir')
    mkdir(childName);
end
cd(childName);
parentName = childName;

childName = sprintf(...
    't-%0.3g_to_%0.3g-step_%0.3g-RandConst_In-mean_%0.3g-stddev_%0.3g',...
    t_init,t_fin,step,in_mean,in_stddev);
childName = fullfile(parentName,childName);
childName = strrep(childName,'.','_');
if ~exist(childName,'dir')
    mkdir(childName);
end
cd(childName);
parentName = childName;

% Create the names for the save files
firingPlotFileName = sprintf(...
    'firingPlot-t_%0.3g_%0.3g-step_%0.3g-E_Percent_%0.2f-LIF_net_%d-RandConstIn-mean_%0.3g-stddev_%0.3g',...
        t_init,t_fin,step,N,E_percent,in_mean,in_stddev);
firingPlotFileName = strrep(firingPlotFileName,'.','_');
firingPlotFileName = strcat(firingPlotFileName,'.jpg');

phasePlotFileName = sprintf(...
    'phasePlot-t_%0.3g_%0.3g-step_%0.3g-E_Percent_%0.2f-LIF_net_%d-RandConstIn-mean_%0.3g-stddev_%0.3g',...
        t_init,t_fin,step,E_percent,N,in_mean,in_stddev);
phasePlotFileName = strrep(phasePlotFileName,'.','_');
phasePlotFileName = strcat(firingPlotFileName,'.avi');
    
% Write the firing plot to an image
saveas(firingFigure,firingPlotFileName);

% and the phase plot to a video.
v = VideoWriter(phasePlotFileName);
open(v);
writeVideo(v,vid);
close(v);

cd(oldFolder);

%##########################################################################