function [ P ] = legendre( x,l )
% legendre: computes the lth legendre polynomial at x
%
% input:
%
% output:
%
%##########################################################################
% Pseudo Code:
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Main Algorithm:
%   ####
%##########################################################################
% Input Format Check:

% Make sure the polynomial order and the evaluation point are given.
if nargin<2,error('Need the polynomial number and the evaluation point.');end

%==========================================================================
% Variable Declarations:

% Precompute the variables needed for the computation
xm = (x-1);
xp = (x+1);

% Initialize the output value
P = 0;

%==========================================================================
% Main Algorithm:

% Run the recursive formula to compute the value of the Legendre polynomial
for i=0:l
    P = P + ((nchoosek(l,i).^(2)) .*((xm).^(l-i)) .*((xp).^(i)));
end
P = P./(2.^(l));
    
%##########################################################################
end

