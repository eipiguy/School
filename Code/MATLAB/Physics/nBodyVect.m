function [ dXV1 ] = nBodyVect( t, XV1, dim, m, G )
%% nBody:
%   [ dx ] = nBody( t, x )
%       Description
% input:
%   t = independent time variable
%   xv = combined position and velocity vector
%     = [ x1; v1; x2; v2; ... ; xn; vn ]
%   m = mass vector
%   G = gravitational constant
% output:
%   dxv = resulting derivative vector
%       = [ dx1; dv1; dx2; dv2; ... ; dxn; dvn ]
%##########################################################################
%% Pseudo-Code:
%   ####
%   Parameters:
%   ====
%   Main ODE:
%   ####
%##########################################################################
%% Parameters:

XV1 = XV1';

% Numer of Masses
n = length(m);

XV = zeros(2*n,dim);

for i=0:n-1
    for j=1:dim
        XV((2*i)+1,j) = XV1((2*dim*i)+j);
    end
end

sysDim = size(XV);

% Radius Matrix:
% distance between the ith and jth points
% and corresponding absolute values
sepR = zeros(n,n,sysDim(2));
R = zeros(n,n);

% Intermediate dv (acceleration) vector
dv = zeros(n,sysDim(2));

% Output Derivative Vector
dXV = zeros(2*n,sysDim(2));
dXV1 = zeros(size(XV1));

%==========================================================================
%% Main ODE:

% Set up the radius matrix
for i=1:n
    for j=i:n
        sepR(i,j,:) = XV((2*j)-1,:) - XV((2*i)-1,:); 
        R(i,j) = sqrt( sum( ( sepR(i,j,:).^(2) ), 3 ) );
        
        sepR(j,i,:) = sepR(i,j,:);
        R(j,i) = R(i,j);
    end
end

% Set up the coresponding absolute value matrix
aR3 = ( abs(R) ).^(-3);
for i=1:n
    for j=1:n
        if ~isfinite(aR3(i,j))
            aR3(i,j) = 0;
        end
    end
end
%display(aR3);

% Calculate the acceleration vector
for i=1:sysDim(2)
    dv(:,i) = -G * ( ( sepR(:,:,i) .* aR3  ) * m );
end

%display(dv);
%display(dXV);
%display(XV);

% Inteleave the result with the velocity for the output
for i=0:(n-1)
    dXV((2*i)+1,:) = XV(2*(i+1),:);
    dXV(2*(i+1),:) = dv(i+1,:);
end

for i=0:n-1
    for j=1:dim
        dXV1((2*dim*i)+j) = dXV((2*i)+1,j);
    end
end

%dXV1 = dXV1';
%display(dXV1);

end    