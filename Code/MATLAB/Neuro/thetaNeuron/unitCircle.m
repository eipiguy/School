function [ ] = unitCircle( h )
%% unitCircle: draws a unit circle with given angle step size between points
%   [ img ] = unitCircle( h )
%       Description
% input:
%   h = angle step size
% output:
%   img = output image
%% Pseudo-Code:
%   ####
%   Main Algorithm:
%   ####
%##########################################################################
%% Main Algortithm:

    theta= 0:h:(2*pi);
    plot(cos(theta),sin(theta));

%##########################################################################
end

