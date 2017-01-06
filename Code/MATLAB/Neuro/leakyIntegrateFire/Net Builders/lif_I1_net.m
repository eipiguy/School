function [ gl, G, El, E, C, ...
    delay, rise, decay, ...
    thresholdMask, threshBases, threshJumps, threshDecays, ...
    refractoryValues, refractoryTimes ] = lif_I1_net( N, pCon )
%% I1_net:  for a network of single population inhibitory neurons 
%   Uses leaky integrate and fire, row to column,
%   to model a single, disconnected neuron,
%   with no input,
%   write the solution to a csv file as
%   [ time (s), transmembrane potential (mV) ],
%   with each entry in a separate row and ordered by time,
%   and then display the results.
%##########################################################################
%% Input:
%   N = number of neurons in population
%
%==========================================================================
%% Output
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
%   Input Check:
%   ====
%   Parameter Assignment:
%   ====
%   Build Connections:
%   ####
%##########################################################################
%% Input Check:

%==========================================================================
%% Parameter Assingment

gl = (0.25e-7)*ones(N,1);   % inhibitory membrane leak conductivity (S)
El = (-68.6)*ones(N,1);     % transmembrane potential (mV) at equilibrium

E = (-80)*ones(N,1);        % inhibitory synaptic equilibrium voltage (mV)
C = (5e-7)*ones(N,1);       % inhibitory membrane capacitance (mF)
%--------------------------------------------------------------------------
thresh = -55*ones(6,1);     % threshold potential (mV)
holdV = -68.6*ones(6,1);    % firing refractory hold potential (mV)
holdT = 1.5e-3*ones(6,1);   % time of post fire refractory period (s)
%--------------------------------------------------------------------------
delay = 0.2;
rise = 0.2;
decay = 0.75;

%##########################################################################
end

