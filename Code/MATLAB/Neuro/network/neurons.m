function [ El, Er, C, GMax, W ] = neurons( n, ePercent, fsRatio, CP, iCP )
%% synapses:
%   [ E, C, S ] = network( n, ePercent, fsRatio )
%       Description
% input:
%   n = number of neurons in network
%   ePercent = percentage of neurons that are excitatory
%   fsRatio = fast to slow ratio of remaining inhibitory neurons
% output:
%   L =
%   El =
%   Er = 
%   C = 
%   G = 
%##########################################################################
%% Pseudo-Code:
%   ####
%   Parameters:
%   ====
%   Output Formatting:
%   ####
%##########################################################################
%% Parameters:

lRest = -68;
elCond = 0.000000025;
ilCond = 0.00000002;

eRev = 0;
eCapa = 0.0000005;

iRev = -80;
iCapa = 0.0000002;

eeCond = 0.000003;
eiCond = 0.000005;
ieCond = 0.000004;
iiCond = 0.000001;

El = zeros(n,1);
Er = zeros(n,1);
C = zeros(n,1);
GMax = zeros(n,n);
W = zeros(n,n);

%==========================================================================
%% Main ODE:

for i=1:n
    
    % Set all the receivers for neuron i
    for j=1:n
        GMax(i,j) = 0;
    end
    
    El(i) = lRest;
    
    % Set the physical characteristics of excitatory vs. inhibitory
    
    % Excitatory indeces
    if (i/n) <= ePercent
        GMax(i,i) = elCond;
        Er(i) = eRev;
        C(i) = eCapa;
        
        for j=1:n
            if (j/n) <= ePercent
                GMax(i,j) = eeCond;
                if rand <= CP
                    W(i,j) = 1;
                end
            else
                GMax(i,j) = eiCond;
                if rand <= iCP
                    W(i,j) = 1;
                end
            end
        end

    % Inhibitory indeces
    else
        Er(i) = iRev;
        C(i) = iCapa;
        
        for j=1:n
            if (j/n) <= ePercent
                GMax(i,j) = ieCond;
                if rand <= CP
                    W(i,j) = 1;
                end
            end
            if (j/n) <= ePercent
                GMax(i,j) = iiCond;
                if rand <= iCP
                    W(i,j) = 1;
                end
            end 
        end
    end
end

%##########################################################################
end