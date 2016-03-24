function [ U ] = cholesky( A )
% cholesky: Cholesky decomposition of a symmetric matrix
%   [ U ] = cholesky( A )
%       uses Cholesky decomposition
%       to decompose A = U^{T}U
%           U = Upper triangular
% input:
%   A = symmetric input matrix
% output:
%   U = upper triangular factor of A
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

% make sure an input matrix and product solution were given
if nargin < 1, error('Please input a matrix.');end

% set the size of the matrix 
% and check to make sure the matrix is square
[m,n] = size(A);
if m~=n, error('Input matrix is not square.');end

% check to make sure the matrix is symmetric
if ~issymmetric(A), error('Input matrix is not symmetric');end
%==========================================================================
% Variable Declarations:

%==========================================================================
% Main Algorithm:

%##########################################################################
end

