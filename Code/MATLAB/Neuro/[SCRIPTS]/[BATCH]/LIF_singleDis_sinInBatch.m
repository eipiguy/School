%% LIF_singleDis_sinIn:
%   Uses the LIF neuron model 
%   for a single, disconnected neuron
%   with a sinusoidal input,
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
%   aIn = amplitute of input sine wave 
%           in miliAmps per centimeter^2 (mA/cm2)
%   fIn = angular frequecy of input sine wave in radians per second (rad/s)
%   pIn = phase of input sine wave in radians (rad)
%   ----
%   In(t) = input control current across membrane 
%           in milliAmps per square centimeter (mA/cm2)
%   ----
%   g = membrane conductance in miliSiemens per centimeter^2 (mS/cm2)
%   E = equilibrium transmembrane potential in miliVolts (mV)
%   C = membrane capacitence in miliFarads per centimeter^2 (mF/cm2)
%   ----
%   thresh = threshold of transmembrane potential in miliVolts (mV)
%               that produces firing behavior
%   holdV = transmembrane potential in miliVolts (mV) 
%               after refractory period
%   holdT = time of firing refractory period in seconds(s)
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
step = 1e-2;    % time interval step size in ``delta" seconds (ds)
%--------------------------------------------------------------------------
V_Init = -60;   % transmembrane potential (mV) at time t_init
%--------------------------------------------------------------------------
aMin = 10;  % amplitude of current into the membrane (mA/cm^2)
aMax = 50;
aStep = 10;
fMin = 0;
fMax = 30;   % current frequency into the membrane (rad/s)
fStep = 10;
pIn = 0;    % phase of current into membrane (rad)
%--------------------------------------------------------------------------
g = 0.25;   % membrane conductivity (mS/cm2)
E = -68.6;  % transmembrane potential (mV) at equilibrium
C = 0.1;    % membrane capacitance (mF/cm2)
%--------------------------------------------------------------------------
thresh = -55;   % threshold potential (mV)
holdV = E;      % firing refractory hold potential (mV)
holdT = 1.5e-3; % time of post fire refractory period (s)
%--------------------------------------------------------------------------
delay = 1.5e-1; % delay before firing triggers synapse reaction
rise = 0.2;     % exponential riseing time constant for the gate function
decay = 0.5;    % exponential decay time constant for the gate function
%==========================================================================
for batchA = aMin:aStep:aMax
    for batchF = fMin:fStep:fMax
%==========================================================================
%% Setup:

        In =@(t) sinIn(t,batchA,batchF,pIn);  % input control current (mA/cm2), 
                                        % dependent on time
        
        dVS =@(t,VS,fireT) lif( t,VS,fireT,delay,rise,decay,In,g,E,C );

%--------------------------------------------------------------------------
%% Solve:

        [ tOut, VOut, fireMask, threshMask ] = eul_resThresh( ...
                                dVS, [ t_init, t_fin ], [ V_Init, 0 ], ...
                                step, [ 1, 0 ], thresh, holdV, holdT );

%--------------------------------------------------------------------------
%% Display:

        m = length(V_Init);
        firingTimes = nan(length(tOut),m);
        for j=1:m
            subPass = tOut;
            subPass(~fireMask(:,j)) = nan;
            firingTimes(:,j) = subPass';
        end

%display(firingTimes);

        clf;
        cla;
        
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
            'firingPlot-singleDisconnectedLIF-sinusoidalInput_a%d_f%d.jpg',...
            batchA,batchF);
        phasePlotFileName = sprintf(...
            'phasePlot-singleDisconnectedLIF-sinusoidalInput_a%d_f%d.avi',...
            batchA,batchF);
        
        saveas(gcf,firingPlotFileName);

        v = VideoWriter(phasePlotFileName);
        open(v);
        writeVideo(v,vid);
        close(v);

    end
end
                                        
%##########################################################################