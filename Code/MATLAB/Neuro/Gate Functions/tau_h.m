function [ tau ] = tau_h( a, b )
%% h_inf: gating function based on other gating functions alpha and beta
%   Detailed explanation goes here
%##########################################################################
%% Input:
%   a = alpha gating variable
%   b = beta gating variable
%==========================================================================
%% Output:
%   tau = value of gating variable
%##########################################################################
%% Pseudocode:
%   ####
%   Main Equation:
%   ####
%##########################################################################
%% Main Equation:

tau = (a+b).^(-1);

%##########################################################################
end

