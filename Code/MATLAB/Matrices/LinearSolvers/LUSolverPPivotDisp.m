function [ X, P, L, U ] = LUSolverPPivotDisp( A, B )
% LUSolverPPivotDisp: LU decomposition solver with partial pivoting
%   [ L, U, P ] = LUPPivotDisp( A, B )
%       uses Gauss elimination with partial pivoting
%       to pivot and decompose PA = LU
%           P = Pivot matrix
%           L = lower triangular
%           U = Upper triangular
%       and then solve the linear system AX=B
%       displaying results as you go
% input:
%   A = input matrix
%   B = set of right hand side vectors 
% output:
%   X = Solutions to the different systems of equations
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
%       Forward Elimination: (2/3)n^3 -(5/6)n^2 +(1/6)n
%           ****
%           Partial Pivoting:
%           ****
%       ----
%       Seperation:
%       ====
%       Forward Substitution: (nB)[(3/2)n^{2}-(1/2)n-1]
%       ----
%       Backward Substitution: (nB)[(1/2)n^{2}+(3/2)n-2]
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

P = eye(n); % initial pivoting matrix before starting algorithm

%==========================================================================
% Main Algorithm:

display(A);   % display the starting matrix

%--------------------------------------------------------------------------
% Forward Elimination:

% for each pivot column
for k = 1:(n-1)
    
    %**********************************************************************
    % Partial Pivoting:
    
    % find the maximum leading term in the current pivot column
    [~,i] = max(abs(A(k:n,k)));
    % if it is not the current row
    if i ~= 1
        % switch the current row and the row with the largest leading term 
        A([k,i+k-1],:)=A([i+k-1,k],:);
        B([k,i+k-1],:)=B([i+k-1,k],:);
        P([k,i+k-1],:)=P([i+k-1,k],:);

        display([P,A]);   % display the result of the pivot
    end    
    %**********************************************************************
    
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
%
% Total Flops: sum_{i=2}^{n}[ (nB)[2(i-1)-1] +(nB) ] 
%   = (nB)n(n-1) = (nB)[ n^2 -n ]

D = zeros(size(B));     % Set up our intermediary solution
D(1,:) = B(1,:);        % Find the initial values to start substitution

for i=2:n
    D(i,:) = B(i,:)-(L(i,1:(i-1))*D(1:(i-1),:));
end
% (n-1) cycles in i,
%   (nB)[2(i-1)-1] flops for producing the variables to subtract:
%       (nB) columns in B, each produced with:
%           (i-1) multiplications and 
%           (i-1)-1 additions
%   (nB) flops for subtracting the row to find the next iteration's values

display(D);

%--------------------------------------------------------------------------
% Back Substitution:
%
% Total Flops: sum_{1}^{n-1}[ (nB)[2(n-i)-1] +(nB) +(nB) ]
% = (nB)(n-1)(n+1) = (nB)(n^2-1)

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
%   (nB)[2(n-i)-1] flops for preparing the previous variables to subtract
%       (nB) columns in B, each produced with:
%           (n-i) multiplications and 
%           (n-i)-1 additions
%   (nB) flops for subtracting them from the intermediate solution
%   (nB) flops for scaling the results

%##########################################################################