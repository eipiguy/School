function [ dx ] = nBody( t, x, m, G )
%% dualBody:
%   [ dx ] = nBody( t, x )
%       Description
% input:
%   t = independent time variable
%   x = combined position and velocity vector
%     = [ x1; v1; x2; v2; ... ; xn; vn ]
%   m = mass vector
%   G = gravitational constant
% output:
%   dx = 
%##########################################################################
%% Pseudo-Code:
%   ####
%   Parameters:
%   ====
%   Main ODE:
%   ####
%##########################################################################
%% Parameters:

n = length(m);

r = zeros(n);
for i=1:n
    r(:,i) = x() - x(i);
end
%==========================================================================
%% Main ODE:

dx = zeros(2*n);

end



    