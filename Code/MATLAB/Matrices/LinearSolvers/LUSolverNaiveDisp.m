function [ L, U, x ] = LUSolverNaiveDisp( A, b )
% LUSolverNaiveDisp: LU decomposition solver
%   [ L, U, P ] = LUSolverNaiveDisp( A, B )
%       uses naive Gauss elimination
%       to pivot and decompose PA = LU
%           P = Pivot matrix
%           L = lower triangular
%           U = Upper triangular
%       and then solve the linear system 
%       displaying results as you go
% input:
%   A = input matrix
%   B = set of right hand side vectors 
% output:
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
%       ----
%       Final Seperation:
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
        A(i,k+1:n) = A(i,(k+1):n) - factor*A(k,(k+1):n);
    end
    % (n-k) cycles in i,
    % 1 flop per for finding elimination row factor, 
    % 2(n-k) flops per with:
    %   n-k for multiplying factor through pivot row and
    %   n-k for subtracting the resulting row from the elimination row
    % Total Flops:(n-k)[2(n-k)+1]
    
    display(A);   % display after each pivot column elimination
end
% (n-1) cycles in k
% Total Flops: sum_{k=1}^{n-1} (n-k)[2(n-k)+1] = 
%   = n(n-1)(4n+5)/6 = (2/3)n^3 +(1/6)n^2 -(5/6)n


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
% (n-1) cycles in i
% No flops, only reassignment.

U = A;

%==========================================================================
% Forward Substitution:

d = zeros(n,1);
d(1) = b(1);
for i=2:n
    d(i) = b(i)-(L(i,1:(i-1))*d(1:(i-1)));
end

display(d);

%--------------------------------------------------------------------------
% Back Substitution:

x = zeros(n,1);
x(n) = d(n)/U(n,n);
% moving backwards through the solution components
for i = n-1:-1:1
    % solve for the next component up
    % that depends only on previously computed components
    x(i) = (d(i) - U(i,(i+1):n)*x((i+1):n))/U(i,i);
end


%##########################################################################