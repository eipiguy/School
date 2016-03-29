function [ x ] = gaussSeidelRelaxDisp( A, b, l, es, maxit )
% gaussSeidelRelaxDisp: Gauss-Seidel naive iterative solver
%   [ x ] = gaussSeidelRelaxDisp( A, b, l, es, maxit )
%       uses the Gauss-Seidel iterative method with relaxation
%       to solve the linear system Ax=b
%       by recursively solving for each component in x
%       adjusting the results using relaxation
%       and displaying results as you go
% input:
%   A = input matrix
%   b = right hand side vector
%   l = relaxation constant: 
%       0-1:underrrelaxed, 1-2: over relaxed, (default = 1 = unchanged)
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
if nargin<2, error('coefficients and right hand side required'); end

if isempty(l) || l==0 || l==1
    error('use naive method without relaxation');
end

% Set defaults for maximum iterations and relative error
% if they aren't explicitly given
if nargin<5 || isempty(maxit), maxit=50; end
if nargin<4 || isempty(es), es=0.00001; end

% Check the size of A to make sure it's square
[m,n] = size(A);
if m~=n, error('Matrix A must be square'); end

%==========================================================================
% Variable Declarations:

% Set up the main matrix for calculation coefficients (size only).
C = A;
d = b;

% Set up the initial solution and relative errors
x = zeros(n,1);
ea = ones(n,1);

% divide every row in the matrix by the corresponding diagonal entry of A
for i = 1:n
    C(i,1:n) = C(i,1:n)/A(i,i);
end

% divide the diagonal entries through the rhs as well
for i = 1:n
    d(i) = b(i)/A(i,i);
end

% Clear the diagonal of the calculation matrix.
for i = 1:n
    C(i,i) = 0;
end

display(d); % display the adjusted rhs vector once it has been scaled
display(C); % display the calculation matrix once it is finished

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
            x(i) = l*x(i)+(1-l)*xold(i);
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

