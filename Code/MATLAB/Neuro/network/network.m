function [ dV ] = network( t, V, E, C, S, In )
%% network:
%   [ dV ] = network( t, V, E, C, S, In )
%       Description
% input:
%   t = independent time variable
%   V = neuron voltage vector
%   E = equilibrium potential vector
%   C = capacitence vector
%   S = synapse connection matrix
%   In = input current function vector
% output:
%   dV = potential change vector
%##########################################################################
%% Pseudo-Code:
%   ####
%   Parameters:
%       G = synapse gating matrix
%   Main ODE:
%   ####
%##########################################################################
%% Parameters:

% Synapse gating matrix
G = zeros(size(S));

dV = (((S.*G)*(V-E)) + In(t) )./(-C);
end