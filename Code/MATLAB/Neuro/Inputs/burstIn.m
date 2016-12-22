function [ InV ] = burstIn( t, V, ti, tf )
%% bburstIn: Single time variable in, burst vector out
%   A burst function that sends a zero vector for a given period of time, 
%   and then between times in vectors ti and tf, outputs the corresponding
%   constant value from V.
%##########################################################################
%% Input:
%   t = intput independent variable
%   V = constant output vector
%   ti = burst start time vector
%   tf = burst ending time vector
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

m = length(V);
InV = zeros(m,1);
for i=1:m
    if t < ti(i) || tf(i) < t
        InV(i) = 0;
    else
        InV(i) = V(i);
    end
end
    
%##########################################################################
end