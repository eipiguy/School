function [ InV ] = delayIn( t, V, t0 )
%% delayIn: Single time variable in, delayed constant vector out
%   A delay function that sends a zero vector for a given period of time, 
%   and then after time in vector t0, outputs the corresponding constant
%   value from V.
%##########################################################################
%% Input:
%   t = intput independent variable
%   V = constant output vector
%   t0 = delay start time vector
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
    if t0(i) <= t
        InV(i) = V(i);
    else
        InV(i) = 0;
    end
end
    
%##########################################################################
end