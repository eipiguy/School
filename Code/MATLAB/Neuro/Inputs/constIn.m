function [ InV ] = constIn( t, V )
%% constIn: Single time variable in, constant vector out
%   A constant vector function that sends a a given vector out for every
%   time value.
%##########################################################################
%% Input:
%   t = intput independent variable
%   V = constant output vector
%==========================================================================
%% Output:
%   I = Resulting vector out at a given time
%##########################################################################
%% Pseudocode:
%   ####
%   Main Algorithm:
%   ####
%##########################################################################
%% Main Algorithm:

InV = V*ones(size(t));
    
%##########################################################################
end