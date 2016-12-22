function [ a ] = alpha_h( V )
%% alpha_h: gating function based on transmembrane potential
%   alpha_h gating function used by Hodgkin-Hukley model
%##########################################################################
%% Input:
%   V = transmembrane potential in miliVolts (mV)
%==========================================================================
%% Output:
%   a = value of gating variable
%##########################################################################
%% Pseudocode:
%   ####
%   Main Equation:
%   ####
%##########################################################################
%% Main Equation:

a = 0.07*exp( -V/20 );

%##########################################################################
end