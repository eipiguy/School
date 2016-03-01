function [ f ] = freefallDrag( cd )
% freefallDrag
%
% This function takes a given drag coefficient,
% and given the values for gravitational acceleration,
% the object's mass, and time of free fall,
% returns the velocity of such an object minus the desired value.

g = 9.81;
m = 80;
t = 4;
v = 36;

f = sqrt(g*m/cd).*tanh(sqrt(g*cd./m)*t)-v;

end

