function [ dV ] = network( t, V, El, Er, C, GMax, W, In )
%% network:
%   [ dV ] = network( t, V, L, E, C, S, In )
%       Description
% input:
%   t = independent time variable
%   V = neuron voltage vector
%   El = receiver neurons' equilibrium voltage
%   Er = reversal potentials as senders
%   C = capacitence vector
%   GMax = synapse conductance matrix
%   W = synapse coupling strength matrix
%   In(t) = input current function vector
% output:
%   dV = potential change vector
%##########################################################################
%% Pseudo-Code:
%   ####
%   Parameters:
%   ====
%   Main ODE:
%   ####
%##########################################################################
%% Parameters:

m = length(V);  % number of neurons

EM = eye(m);    % leak conductance matrix
W = ones(m,m);  % coupling strength matrix
S = zeros(m,m);  % synapse gate matrix



% Equilibrium Matrix
for i=1:m
    EM(i,i) = GMax(i,i);
    GMax(i,i) = 0;
end
%==========================================================================
%% Main ODE:

dV = ( In(t) - ( GMax*(V'-Er) ) - ( EM*(V'-El) ) )./C;

end