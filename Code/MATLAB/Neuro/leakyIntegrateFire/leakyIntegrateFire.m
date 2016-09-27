function [ dV ] = leakyIntegrateFire( t, V, ePercent, fSRatio, In )
%% leakyIntegrateFire: LIF population ODE system
%   [ dtheta ] = thetaNeuron( t, theta, I )
%       Description
% input:
%   t = independent time variable
%   V = current voltage for each neuron:
%       excitatory neurons have the first indeces,
%       then fast inhibitory, then slow inhibitory
%   ePercent = percentage of the neuron population that is excitatory
%   fSRatio = the number of fast neurons per slow in the population
%   I = applied current for each neuron
% output:
%   dtheta = change in angle on the unit circle
%##########################################################################
%% Pseudo-Code:
%   ####
%   Parameters:
%   ====
%   Neuron Population:
%   ====
%   ODE Equation:
%   ####
%##########################################################################
%% Parameters:

% Rest Potentials (mV)
eL = -68;   % leak
eE = 0;     % excitatory
eI = -80;   % inhibitory

% Capacitence (mF)
cE = 5*(10^(-7));   % excitatory
cL = 5*(10^(-7));   % inhibitory

% Leak Conductance (S)
glE = 0.25*(10^(-7));   % exciatory
glE = 0.25*(10^(-7));   % inhibitory

% Refreactory Period (ms)
rE = 3;     % excitatory
rI = 1.5;   % inhibitory

% Threshold Refractory Effects
tJ = 10;    % threshold jump (mv)
tTC = 10;   % threshold decreasing time constant (ms)

%==========================================================================
%% Neuron Population:

nPop = length(V);   % the number of neurons in the population

% Set the inputs for the respective neuron types
for i=1:nPop
    
    % First set the extitatory population
    if (i/nPop) <= ePercent
        
        % Set the leak influence and 
        % the direct input's influence
        dV(i) = (gIn.*In(i)) - (gL.*(V(i)-eL));
        
        % Set any excitatory influence
        if true
            dV(i) = dV(i) - (gE);
        end
        
        % Set any inhibitory influence
        if true
            
        end       
        
    % Then set the fast inhibitory population
    else
        
    end
end

%==========================================================================
%% ODE Equation:

    dV = -V + In(t);        

%##########################################################################
end