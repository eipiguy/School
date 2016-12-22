function [ I ] = sinIn( t, a, f, p )
%% sinIn: Single time variable in, sin vector of given components out
%   Outputs a time dependent set of sine waves with the given 
%   amplitutes, angular frequencies, and phases.
%##########################################################################
%% Input:
%   t = intput independent variable
%   a = amplitute vector
%   f = angular frequency vector
%   p = phase vector
%==========================================================================
%% Output:
%   I = Output sine value at time t
%##########################################################################
%% Pseudocode:
%   ####
%   Main Algorithm:
%   ####
%##########################################################################
%% Main Algorithm:

    % Check the number of variables
    % and initialize the output vector
    m = length(a);
    I = zeros(m,1);
    
    % Calculate the sine value for each component and add it to the output
    for i=1:m
        I(i) = a(i)*sin( (f(i)*t) + p(i) );
    end

%##########################################################################
end