function x = gaussPPivot( A, b )
% gaussPPivot: linear solver for gaussian elimination with partial pivoting
%   [ x ] = gaussPPivot( A, b )
%       uses Gauss elimination with partial pivoting 
%       to solve the linear system Ax = b
% input:
%   A = linear system coefficients
%   b = product vector
% output:
%   x = solution vector to produce product with given coefficient matrix
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
%       Back Substitution:
%   ####
%##########################################################################
% Input Format Check:

% make sure an input matrix and product solution were given
if nargin < 2, error('Not enough inputs.');end

% set the size of the matrix 
% and check to make sure the matrix is square
[m,n] = size(A);
if m~=n, error('Input matrix is not square.');end

%==========================================================================
% Variable Declarations:

aug = [A b];    % the augmented matrix for elimination to solve
nb = n+1;       % the number of columns in the augmented matrix

%==========================================================================
% Main Algorithm:

%--------------------------------------------------------------------------
% Forward Elimination:

% for each pivot column
for k = 1:n-1
    
    %**********************************************************************
    % Partial Pivoting:
    
    % find the maximum leading term in the current pivot column
    [~,i] = max(abs(aug(k:n,k)));
    % if it is not the current row
    if i ~= 1
        % switch the current row and the row with the largest leading term 
        aug([k,i+k-1],:)=aug([i+k-1,k],:);
    end
    
    %**********************************************************************
    
    % for each row past the diagonal 
    for i = k+1:n
        
        % eliminate the leading term using row operations:
        
        % find the factor that changes the column's leading term
        % to match the leading term of the current elimination row
        factor = aug(i,k)/aug(k,k);
        
        % multiply the factor through the pivot row
        % and subtract the result from the current elimination row
        aug(i,k:nb) = aug(i,k:nb) - factor*aug(k,k:nb);
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
% Back Substitution:

% set up the place to store the final solution
x = zeros(n,1);

% set the last variable in the solution
% since it is given explicitly after elimination
x(n) = aug(n,nb)/aug(n,n);

% moving backwards through the solution components
for i = n-1:-1:1
    % solve for the next component up
    % that depends only on previously computed components
    x(i) = (aug(i,nb) - aug(i,i+1:n)*x(i+1:n))/aug(i,i);
end

%##########################################################################