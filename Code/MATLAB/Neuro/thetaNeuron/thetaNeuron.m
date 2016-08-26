function [ dtheta ] = thetaNeuron( t, theta, I )
%% thetaNeuron: single parameter ODE for neuron action potential
%   [ dtheta ] = thetaNeuron( t, theta, I )
%       Basic ODE for the theta-neuron (Ermentrout-Kopell) model.
%       Plots action potential as y-component on the unit circle;
%       theta = pi is considered a neuron burst, while the input
%       is in terms of current.
% input:
%   t = independent time variable
%   theta = dependent angle on the unit circle
%       (corresponds to action potential)
%   I = input function (applied current)
% output:
%   dtheta = change in angle on the unit circle
%##########################################################################
%% Pseudo-Code:
%   ####
%   ODE Equation:
%   ####
%##########################################################################
%% ODE Equation:

    dtheta = 1 - cos(theta) + (( 1 + cos(theta) ).*(I(t)));        

%##########################################################################
end