function [ dVS ] = lif_singleInhibit( t, VS, lastFire, delay, rise, decay, In, g, El, Ei, C )
%% lif: leaky integrate and fire ODE, single neuron, row in/row out
%   [ dVS ] = lif( t, VS )
%       Uses a column vector representing transmembrane potential and 
%       a gating variable to model the activity of a single neuron.
%##########################################################################
%% Variable Declaration:
% input:
%   t = independent column vector time in seconds
%   VS = dependent row vector
%       [   V = transmembrane potential in miliVolts (mV),
%           S = gating variable   ]
%   ----
%   In(t) = Input control current function,
%           current in microAmps/cm^(2) (mA/cm2),
%           at a given time t in seconds (s)
%   ----
%   lastFire = firing times in seconds (s)
%   delay = gating reaction delay in seconds(s)
%   rise = constant that determines gate opening rate
%   decay = constant that determines gate closing rate
%   ----
%   g = membrane conductance in miliSiemens/centimeter^(2) (mS/cm2)
%   gi = inhibitor synapse conductance in miliSiemens/cm^(2) (mS/cm2)
%   El = leak equilibrium potential in miliVolts (mV)
%   Ei = inhibitor synapse equilibrium potential in miliVolts (mV)
%   C = membrane capacitane in miliFarads/centimeter^(2) (mF/cm2)
% output:
%   dVS = vector that indicates direction of change in phase space
%##########################################################################
%% Pseudocode:
%   ####
%   Main Equation:
%   ####
%##########################################################################
%% Variables:

% Record the number of dimensions
m = size(lastFire,2);

%==========================================================================
%% Transmembrane Potential:

% Set the voltage change
dVS(1) = ( In(t) - (g*(VS(1)-El)) - (gi*(VS(1)-Ei)) ) /C;

%--------------------------------------------------------------------------
%% Synaptic Gate:

% Set the decay term of the gating variable
dVS(2) = -VS(2)./decay;

% Check the firing times to see if any activate their spike gate
times = lastFire(:,1);
times = times(~isnan(times));
tL = length(times);

spikeFlag = false;
if ~isempty(times)
    for i=1:tL
        if (( 0 < (t-times(i)-delay) ) && ( (t-times(i)-delay) < 1 ))
            spikeFlag = true;
        end
    end
end

if spikeFlag, dVS(2) = dVS(2) + ( (1-VS(2))./rise ); end

%##########################################################################
end