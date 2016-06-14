function [ config, acceptRate ] = nsphBox( dim, box, num, rad )
%% nsphBox: equiprobable configuration generator
%   [ config, acceptRate ] = nsphBox( dim, box, no, rad )
%       Randomly generates a set of num coordinate touples
%       placed such that they are at least 2*rad apart pairwise,
%       and at least rad away from the edges of the given boundaries.
%       Uses a tablua rasa rule to ensure equiprobability of every state,
%       and returns the acceptance rate given random placements of points.
% input:
%
% output:
%
%% Pseudocode:
% Input Format Check:
%% Input Format Check:
