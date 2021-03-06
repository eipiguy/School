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

[m,n]=size(A);
U = zeros(m,n);

%==========================================================================
% Main Algorithm:

U(1,1)=sqrt(A(1,1));
for j=2:n
    U(1,j) = A(1,j)/U(1,1);
end

for i=2:n
    U(i,i) = sqrt( A(i,i) - (U(1:(i-1),i)'*U(1:(i-1),i)) );
    for j=(i+1):n
        U(i,j) = sqrt(A(i,j) - (U(1:(i-1),i)'*U(1:(i-1),j)));
        display([i,j]);
        display([U(1:(i-1),i)]);
    end
end
%##########################################################################
end

