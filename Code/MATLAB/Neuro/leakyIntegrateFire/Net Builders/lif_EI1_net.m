function [ gl, G, El, E, C, ...
    delay, rise, decay, ...
    thresholdMask, threshBases, threshJumps, threshDecayTConsts, ...
    refractoryValues, refractoryTimes ] = lif_EI1_net( N, E_percent, ...
                                            EIcon_percent_matrix )
%% EI_net:  for a network of single population inhibitory neurons 
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
%   EI_ratio = ratio of excitatory to inhibitory neurons
%   IFS_ratio = ratio of fast to slow inhibitory neurons
%==========================================================================
%% Output
%   gl = [ gl_1; gl_2; ...; gl_N ]
%   gl_i = neuron i's membrane leak conductance 
%           in miliSiemens per centimeter^2 (mS/cm2)
%   G = [   G_1_1,  ...,    G_1_N;
%           ...,    ...,    ...;
%           G_N_1,  ...,    G_N_N   ]
%   G_i_j = maximum weighted net conductance of the synapse 
%           out of neuron i and in to neuron j.
%   El = [ El_1; El_2; ...; El_N ]
%   El_i = neuron i's leak equilibrium potential in miliVolts (mV)
%   E = [  ]
%   E = 
%   C = [ C1; C2; ...; Cn ]
%   Ci = membrane capacitence in miliFarads per centimeter^2 (mF/cm2)
% -------------------------------------------------------------------------
%   thresholdMask = [ thresh1; thresh2; thresh3; ...; threshn ]
%   thresholdMask = threshold of transmembrane potential in miliVolts (mV)
%               that produces firing behavior
%   threshBases = [ holdV1; holdV2; ...; holdVn ]
%   threshBases = transmembrane potential in miliVolts (mV) 
%               after refractory period
%   threshJumps = [  ]
%   threshJumps = 
%   threshDecays = [  ]
%   threshDecays =
%--------------------------------------------------------------------------
%   refractoryValues = [  ]
%   refractoryValues =
%   refractoryTimes = [ holdT1; holdT2; ...; holdTn ]
%   refractoryTimes = time of firing refractory period in seconds(s)
%--------------------------------------------------------------------------
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
for i=1:N
    if i >= (N*E_percent)
        gl(i) = 0.2e-7;
    end
end

G = zeros(N,N);
for i=1:N
    for j=1:N
        if i <= (N*E_percent)
            if j <= (N*E_percent)
                if rand() < EIcon_percent_matrix(1,1) && i~=j;
                    G(i,j) = 3e-6;
                end
            else
                if rand() < EIcon_percent_matrix(1,2) && i~=j;
                    G(i,j) = 4e-6;
                end
            end
        else
            if j <= (N*E_percent)
                if rand() < EIcon_percent_matrix(2,1) && i~=j;
                    G(i,j) = 0.5e-6;
                end
            else
                if rand() < EIcon_percent_matrix(2,2) && i~=j;
                    G(i,j) = 1e-6;
                end
            end
        end
    end
end

El = (-68.6)*ones(N,1);     % transmembrane potential (mV) at equilibrium

E = zeros(N,1);
for i=1:N
    if i > (N*E_percent)
        E(i) = -80;
    end
end

C = (0.25e-7)*ones(N,1);       % inhibitory membrane capacitance (mF)
for i=1:N
    if i > (N*E_percent)
        C(i) = 0.2e-7;
    end
end
%--------------------------------------------------------------------------
thresholdMask = zeros(1,2*N);
for i=1:N
    thresholdMask((2*i)-1) = true;
end

threshBases = -55*ones(1,2*N);            % threshold potential (mV)
threshJumps = 10*ones(1,2*N);           
threshDecayTConsts = 10*ones(1,2*N);    
%--------------------------------------------------------------------------
refractoryValues = nan(2*N,1);    % firing refractory hold potential (mV)
for i=1:N
    refractoryValues((2*i)-1) = El(i);
end

refractoryTimes = 3*ones(2*N,1);   % time of post fire refractory period (s)
for i=1:N
    if i >= (N*E_percent)
        refractoryTimes((2*i)-1) = 1.5;
    end
end
%--------------------------------------------------------------------------
delay = (0.5)*ones(N,1);
rise = (0.2)*ones(N,1);
for i=1:N
    if i >= (N*E_percent)
        rise(i) = 0.6;
    end
end
decay = 10*ones(N,1);
%##########################################################################
end

