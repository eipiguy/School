function [ X, L, U ] = LUSolverNaive( A, B )
% LUSolverNaiveDisp: naive LU decomposition solver
%   [ X, L, U ] = LUSolverNaive( A, B )
%       uses naive Gauss elimination
%       to decompose PA = LU
%           L = lower triangular
%           U = Upper triangular
%       and then solve the linear system AX=B
% input:
%   A = input matrix
%   B = set of right hand side vectors 
% output:
%   X = Solutions to the different systems of equations
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
%       Forward Elimination: (2/3)n^3 -(5/6)n^2 +(1/6)n
%       ----
%       Seperation:
%       ====
%       Forward Substitution: (nB)((3/2)n^{2}-(1/2)n-1)
%       ----
%       Backward Substitution: (nB)((1/2)n^{2}+(3/2)n-2)
%   ####
%##########################################################################
% Input Format Check:

% make sure an input matrix and product solution were given
if nargin < 2, error('Please input coefficients and rhs solutions.');end

% set the size of the matrix 
% and check to make sure the matrix is square
[m,n] = size(A);
if m~=n, error('Input matrix is not square.');end

%==========================================================================
% Variable Declarations:

%==========================================================================
% Main Algorithm:

display(A);   % display the starting matrix

%--------------------------------------------------------------------------
% Forward Elimination:

% for each pivot column
for k = 1:(n-1)
    % for each row past the diagonal 
    for i = k+1:n
        
        % eliminate the leading term using row operations:
        
        % find the factor that changes the column's leading term
        % to match the leading term of the current elimination row
        factor = A(i,k)/A(k,k);
        
        % set the leading term (that normally is eliminated)
        % instead to its corresponding factor to store until needed
        A(i,k)=factor;
        
        % multiply the factor through the pivot row
        % and subtract the result from the current elimination row
        A(i,(k+1):n) = A(i,(k+1):n) - factor*A(k,(k+1):n);
    end
    % (n-k) cycles in i,
    %   1 flop per for finding elimination row factor, 
    %   2(n-k) flops per with:
    %       n-k for multiplying factor through pivot row and
    %       n-k for subtracting the resulting row from the elimination row
    % Total Flops:(n-k)[2(n-k)+1]
    
    display(A);   % display after each pivot column elimination
end
% (n-1) cycles in k
% Total Flops: sum_{k=1}^{n-1} (n-k)[2(n-k)+1] = 
% = n(n-1)(4n-1)/6 = (2/3)n^3 -(5/6)n^2 +(1/6)n


%--------------------------------------------------------------------------
% Seperation:

% set up the place to store L
L = eye(n);

% copy each column of A below the diagonal into L
% and set the value in A to zero
for i = 1:n-1
    % copy each column of A below the diagonal into L
    L((i+1):n,i) = A((i+1):n,i);
    A((i+1):n,i) = 0; % and set original the value in A to zero
end

U = A;

%==========================================================================
% Forward Substitution:

D = zeros(size(B));     % Set up our intermediary solution
D(1,:) = B(1,:);      % Find the initial values to start substitution

for i=2:n
    D(i,:) = B(i,:)-(L(i,1:(i-1))*D(1:(i-1),:));
end
% (n-1) cycles in i,
%   (nb)(i-1) flops for producing the variables to subtract
%   (nb) flops for subtracting the row to find the next iteration's values
% Total Flops: (n-1)((nB)n+(nB))+(nB)sum_{i=2}^{n}(i-1)
% = (nB)(n-1)((3/2)n+1) = (nB)((3/2)n^{2}-(1/2)n-1)

display(D);

%--------------------------------------------------------------------------
% Back Substitution:

X = zeros(size(B));     % Set up the place to put the final solution
X(n,:) = D(n,:)/U(n,n); % Find the last component to begin substitution
% (nB) flops for dividing through

% moving backwards through the solution components
for i = n-1:-1:1
    % solve for the next component up
    % that depends only on previously computed components
    X(i,:) = (D(i,:) - U(i,(i+1):n)*X((i+1):n,:))/U(i,i);
end
% (n-1) cycles in i,
%   (nB)(n-i) flops for preparing the previous variables to subtract
%   (nB) flops for subtracting them from the intermediate solution
%   (nB) flops for scaling the results
% Total Flops: 2(nB)(n-1) + (nB)[ n(n-1)-sum_{1}^{n-1}i ]
% = (nB)(n-1)((1/2)n+2) = (nB)((1/2)n^{2}+(3/2)n-2)

%##########################################################################