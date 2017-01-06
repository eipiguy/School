function [ dVS ] = lif_net( t, VS, fireTimes, In, ...
                                gl, G, El, E, C, delay, rise, decay )
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
%   fireTimes = firing time vector in seconds (s)
%   ----
%   In(t) = Input control current function,
%           current in microAmps/cm^(2) (mA/cm2),
%           at a given time t in seconds (s)
%   ----
%   gl = [ gl_1; gl_2; ...; gl_m ]
%        membrane leak conductance vector 
%        in miliSiemens/centimeter^(2) (mS/cm2)
%   G =  [ g_11, ..., g_1m;
%          ..., ..., ...;
%          g_m1, ..., g_mm ]
%       weighted maximum conductance matrix of synapses
%       from column neuron to row neuron in miliSiemens/centimeter^(2)
%           note that g is not the conductance 
%           of the individual connections, but a weighted conductance
%           representing the total effects of all factors at the boundary,
%               such as: number of synapses, location of connection, etc...
%   El = [ El_1; El_2; ...; El_m]
%        equilibrium potential vector of membrane leak in miliVolts (mV)
%   E = [ El_1; El_2; ...; El_m]
%        equilibrium potential vector for outbound synapses in miliVolts (mV)
%   C = [ C_1; C_2; ...; C_m ]
%       membrane capacitane vector in miliFarads/centimeter^(2) (mF/cm2)
%   ----
%   delay = gating reaction delay in seconds(s)
%   rise = constant that determines gate opening rate
%   decay = constant that determines gate closing rate
%==========================================================================
%% Output:
%   dVS = vector that indicates slope in phase space
%##########################################################################
%% Pseudocode:
%   ####
%   Main Equation:
%   ####
%##########################################################################
%% Variables:

% Record the number of dimensions
m = length(VS)/2;

% Parse the voltage and gate variables
V = zeros(m,1);
S = zeros(m,1);
lastVFire = zeros(m,1);

for i=1:m
    V(i) = VS((2*i)-1);
    S(i) = VS(2*i);
    fires = fireTimes(~isnan(fireTimes(:,(2*i)-1)),(2*i)-1);
    fires = fires(fires <= t);
    if ~isempty(fires)
        lastVFire(i) = fires(end);
    else
        lastVFire(i) = nan;
    end
end

%==========================================================================
%% Transmembrane Potential:

% Set the voltage change
syn = zeros(m,1);
for i=1:m
    for j=1:m
        syn(i) = syn(i) + (G(i,j).*S(j).*(V(i)-E(j)));
    end
end
dV = (In(t) - (gl.*(V-El)) - syn)./C;

%--------------------------------------------------------------------------
%% Synaptic Gate:
if(S(i)==0)
    dS(i)=0;
else
    % Set the decay term of the gating variable
    dS = -S./decay;
end

% Check the firing times to see if any activate their respective gate
for i=1:m
    if ~isnan(lastVFire(i)) && ~isempty(lastVFire(i))
        if (( 0 < (t-lastVFire(i)-delay(i)) ) && ( (t-lastVFire(i)-delay(i)) < 1 ))
            dS(i) = dS(i) + ( (1-S(i))./rise(i) );
        end
    end
end

%==========================================================================
%% Reparse Output:

dVS = zeros(2*m,1);

for i=1:m
    dVS((2*i)-1) = dV(i);
    dVS(2*i) = dS(i);
end

dVS = dVS';

%##########################################################################
end