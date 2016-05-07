function [ x ] = tridiag( e,f,g,r )
%% tridiag: tridiagonal solver banded system
%   [ x ] = tridiag( e,f,g,r )
%       Uses forward elimination and back substitution
%       to solve a set of tridiagonal systems given by band vectors
% input:
%   e = subdiagonal vectors
%   f = diagonal vectors
%   g = superdiagonal vectors
%   r = right hand side vectors
% output:
%   x = solution vectors
%##########################################################################
%% Pseudo Code:
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations:
%   ====
%   Forward Elimination:
%   ====
%   Back Substitution:
%   ####
%##########################################################################
%% Input Format Check:

% Make sure all inputs are given.
if nargin<4,error('All 3 bands and right hand side required.');end

%==========================================================================
%% Variable Declarations:

n = length(f);

%==========================================================================
%% Forward Elimination:

for k=2:n
    factor = e(k,:)./f(k-1,:);
    f(k,:) = f(k,:) -factor.*g(k-1,:);
    r(k,:) = r(k,:) -factor.*r(k-1,:);
end

%==========================================================================
%% Back Substitution:

x(n,:) = r(n,:)./f(n,:);
for k = (n-1):-1:1
    x(k,:) = (r(k,:)-g(k,:).*x(k+1))./f(k,:);
end

%##########################################################################
end