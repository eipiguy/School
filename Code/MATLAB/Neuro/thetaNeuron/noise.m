function [ e ] = noise( t,varargin )
%% noise: noise in ODE based on time
% input:
%   t = time
% output:
%   e = error in ODE differential at time t
%% Pseudo-Code:
%   ####
%   Main Algorithm:
%   ####
%##########################################################################
%% Main Algorithm:

    e = 2*(rand-0.5);
    
%##########################################################################
end

