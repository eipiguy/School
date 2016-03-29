function [ x ] = gaussSeidelNaiveDisp( A, b, es, maxit )
% gaussSeidelNaiveDisp: Gauss-Seidel naive iterative solver
%   [ x ] = gaussSeidelNaiveDisp( A, b, es, maxit )
%       uses the naive Gauss-Seidel iterative method
%       to solve the linear system Ax=b
%       by recusively solving for each component in x
%       displaying results as you go
% input:
%   A = input matrix
%   b = right hand side vector
%   es = stop criterion based on relative error (default = 0.00001%)
%   maxit = maximum allowed iterations (default = 50)
% output:
%   x = solution to the linear system
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

% Check to make sure there are coefficients and a rhs vector.
if nargin<2, error('coefficients and right hand side required');end

% Set defaults for maximum iterations and relative error
% if they aren't explicitly given
if nargin<4 || isempty(maxit), maxit=50; end
if nargin<3 || isempty(es), es=0.00001; end

% Check the size of A to make sure it's square
[m,n] = size(A);
if m~=n, error('Matrix A must be square'); end

%==========================================================================
% Variable Declarations:

% Set up the main matrix for calculation coefficients for size only.
C = A;
d = b;

% Set up the initial solution and relative errors
x = zeros(n,1);
ea = ones(n,1);

% Set the diagonal of calculation matrix to 0.
for i = 1:n
    C(i,i) = 0;
end

% divide every row in the matrix by the corresponding diagonal entry of A
for i = 1:n
    C(i,1:n) = C(i,1:n)/A(i,i);
end
display(C); % display the calculation matrix once it is finished

% divide the diagonal entries through the rhs as well
for i = 1:n
    d(i) = b(i)/A(i,i);
end
display(d); % display the adjusted rhs vector once it has been scaled

iter = 0;   % Initialize the iteration counter

%==========================================================================
% Main Algorithm:

% for each iterative estimate of x
while(1)
    % save the old estimate
    xold = x;
    % calculate the new estimate 
    % while keeping track of the relative error for each component
    for i = 1:n
        x(i) = d(i) - C(i,:)*x;
        if x(i) ~= 0
            ea(i) = abs((x(i)-xold(i))/x(i)) * 100;
        end
    end
    % incrememnt the iteration counter
    iter = iter+1;
    
    display(iter);
    display(x);  % display the current iteration and estimate
    
    % and break if we're within our tolerances
    if max(ea)<=es || iter >= maxit, break, end
end

%##########################################################################
end

