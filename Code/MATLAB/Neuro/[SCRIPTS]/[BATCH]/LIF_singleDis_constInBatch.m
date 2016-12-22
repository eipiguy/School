%% LIF_singleDis_constInBatch:
%   Uses leaky integrate and fire, row to column,
%   to model a single, disconnected neuron
%   with a variety of constant inputs,
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
%   mV_init = transmembrane potential in milliVolts (mV) at time t_init 
%   ----
%   currentIn = current actively added to membrane
%               in miliAmps per cnetimeter^2 (mA/cm2)
%   ----
%   In(t) = input control current across membrane 
%           in milliAmps per square centimeter (mA/cm2)
%   ----
%   g = membrane conductance in miliSiemens per centimeter^2 (mS/cm2)
%   E = transmembrane potential at equilibrium in miliVolts (mV)
%   C = membrane capacitence in miliFarads per centimeter^2 (mF/cm2)
%   ----
%   delay = 
%   rise = 
%   decay = 
%##########################################################################
%% Pseudocode:
%   ####
%   Parameters:
%   ====
%   Setup:
%   ----
%   Solve:
%   ----
%   Save:
%   ----
%   Display:
%   ####
%##########################################################################
%% Parameters:
t_init = 0;     % initial time (s)
t_fin = 3;    % final time (s)
step = 1e-2;    % time interval step size in seconds (s)
%--------------------------------------------------------------------------
V_Init = -60;   % transmembrane potential (mV) at time t_init
%--------------------------------------------------------------------------
minCurrentIn = -10; % constant current into the membrane (mA/cm2)
maxCurrentIn = 10;
currentInStep = 2;
%--------------------------------------------------------------------------
In =@(t) constIn(t,currentIn);  % input control current (mA/cm2),
                                % dependent on time
%--------------------------------------------------------------------------
g = 0.25;    % membrane conductivity (mS/cm2)
E = -68.6;  % transmembrane potential (mV) at equilibrium
C = 0.1;    % membrane capacitance (mF/cm2)
%--------------------------------------------------------------------------
thresh = -55;   % threshold potential (mV)
holdV = E;      % firing refractory hold potential (mV)
holdT = 1.5e-3; % time of post fire refractory period (s)
%--------------------------------------------------------------------------
delay = 0.15;
rise = 0.2;
decay = 0.5;
%==========================================================================
for batch = minCurrentIn:currentInStep:maxCurrentIn
%==========================================================================
%% Setup:
    
    In =@(t) constIn(t,batch);  % input control current (mA/cm2),
                                    % dependent on time
                                
    dVS =@(t,VS,fireT) lif( t,VS,fireT,delay,rise,decay,In,g,E,C );
    
%--------------------------------------------------------------------------
%% Solve:

    [ tOut, VOut, fireMask, threshMask ] = eul_resThresh( ...
                                dVS, [ t_init, t_fin ], [ V_Init, 0 ], ...
                                step, [ 1, 0 ], thresh, holdV, holdT );
                            
%--------------------------------------------------------------------------
%% Display:

    clf;
    cla;
    
    m = length(V_Init);
    firingTimes = nan(length(tOut),m);
    for j=1:m
        subPass = tOut;
        subPass(~fireMask(:,j)) = nan;
        firingTimes(:,j) = subPass';
    end

    %display(firingTimes);

    vid = phasePlot(tOut,VOut,dVS,firingTimes);
    
    figure
    hold on
    for j=1:m
        fire = j*fireMask(:,j);
        ind = find(fire);
        tFire = tOut(ind);
        iFire = fire(ind);
        firePlot = plot(tFire,iFire,'ro');
    end
    xlim manual;
    xlim([t_init t_fin]);

%--------------------------------------------------------------------------
%% Save:

    firingPlotFileName = sprintf(...
        'firingPlot-singleDisconnectedLIF-constantInput_%d.jpg',batch);
    phasePlotFileName = sprintf(...
        'phasePlot-singleDisconnectedLIF-constantInput_%d.avi',batch);

    saveas(gcf,firingPlotFileName);

    v = VideoWriter(firingPlotFileName);
    open(v);
    writeVideo(v,vid);
    close(v);

%##########################################################################
end