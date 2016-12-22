function [ h ] = h_inf( a, b )
%% h_inf: gating function based on other gating functions alpha and beta
%   Detailed explanation goes here
%##########################################################################
%% Input:
%   a = alpha gating variable
%   b = beta gating variable
%==========================================================================
%% Output:
%   h = value of gating variable
%##########################################################################
%% Pseudocode:
%   ####
%   Main Equation:
%   ####
%##########################################################################
%% Main Equation:

h = a./(a+b);

%##########################################################################
end

