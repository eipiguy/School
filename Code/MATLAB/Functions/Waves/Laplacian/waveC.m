function M = waveC( A,B,c,x )
% wave:
%
% input:
%
% output:
%
%##########################################################################
% Pseudo Code:
%   ####
%   Input Format Check:
%   ==== Variable Declarations:
%   ==== Main Algorithm:
%   ####
%##########################################################################
% Input Format Check:

% Make sure the polynomial order and the evaluation point are given.
if nargin<4,error('Need all 3 coefficients and point of evaluation.');end

% Variable Declarations:===================================================


% Main Algorithm:==========================================================

M = (A.*exp(c.*x)) +(B.*exp(-c.*x));
    
%##########################################################################
end

