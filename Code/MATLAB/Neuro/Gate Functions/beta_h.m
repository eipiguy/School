function [ b ] = beta_h( V )
%% beta_h: gating function based on transmembrane potential
%   beta_h gating function used by Hodgkin-Hukley model
%##########################################################################
%% Input:
%   V = transmembrane potential in miliVolts (mV)
%==========================================================================
%% Output:
%   b = value of gating variable
%##########################################################################
%% Pseudocode:
%   ####
%   Main Equation:
%   ####
%##########################################################################
%% Main Equation:

b = ( 1 + exp( (30-V)./10 ) ).^(-1);

%##########################################################################
end

