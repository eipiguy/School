function [ L, U, P ] = LUPPivot( A )
% LUPPivot: LU decomposition using gaussian elim with partial pivoting
%   [ L, U, P ] = LUPPivot( A )
%       uses Gauss elimination with partial pivoting
%       to pivot and decompose PA = LU
%           P = Pivot matrix
%           L = lower triangular
%           U = Upper triangular
% input:
%   A = input matrix
% output:
%   P = Pivot matrix
%   L = lower triangular factor of A
%   U = upper triangular factor of A
%##########################################################################
% Pseudo Code:
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Main Algorithm:
%       ----
%       Forward Elimination:
%           ****
%           Partial Pivoting:
%           ****
%       ----
%       Final Seperation:
%   ####
%##########################################################################
% Input Format Check:

% make sure an input matrix and product solution were given
if nargin < 1, error('Please input a matrix.');end

% set the size of the matrix 
% and check to make sure the matrix is square
[m,n] = size(A);
if m~=n, error('Input matrix is not square.');end

%==========================================================================
% Variable Declarations:

P = eye(n); % initial pivoting matrix before starting algorithm

%==========================================================================
% Main Algorithm:

%--------------------------------------------------------------------------
% Forward Elimination:

% for each pivot column
for k = 1:n-1
    
    %**********************************************************************
    % Partial Pivoting:
    
    % find the maximum leading term in the current pivot column
    [~,i] = max(abs(A(k:n,k)));
    % if it is not the current row
    if i ~= 1
        % switch the current row and the row with the largest leading term 
        A([k,i+k-1],:)=A([i+k-1,k],:);
        % and do the same to the permutation matrix
        P([k,i+k-1],:)=P([i+k-1,k],:);
    end    
    %**********************************************************************
    
    % for each row past the diagonal 
    for i = k+1:n
        
        % eliminate the leading term using row operations:
        
        % find the factor that changes the column's leading term
        % to match the leading term of the current elimination row
        factor = A(i,k)/A(k,k);
        
        % set the leading term (that normally is eliminated)
        % instead to the factor used to eliminate it
        A(i,k)=factor;
        
        % multiply the factor through the pivot row
        % and subtract the result from the current elimination row
        A(i,k+1:n) = A(i,k+1:n) - factor*A(k,k+1:n);
    end
    % (n-k) cycles in i,
    % 1 flop per for finding elimination row factor, 
    % 2(n+2-k) flops per with:
    %   n+2-k for multiplying factor through pivot row and
    %   n+2-k for subtracting the resulting row from the elimination row
    % Total Flops:(n-k)[1+2(n+2-k)] = (n-k)[2(n-k)+5]
end
% (n-1) cycles in k
% Total Flops: sum_{k=1}^{n-1} (n-k)[2(n-k)+5] = 
%   = n(n-1)(4n+13)/6 = 2n^3/3 +3n^2/2 -13n/6


%--------------------------------------------------------------------------
% Final Seperation:

% set up the place to store L
L = eye(n);

% copy each column of A below the diagonal into L
% and set the value in A to zero
for i = 1:n-1
    % copy each column of A below the diagonal into L
    L(i+1:n,i) = A(i+1:n,i);
    A(i+1:n,i) = 0; % and set original the value in A to zero
end
% (n-1) cycles in i
% No flops, only reassignment.

U = A;
%##########################################################################