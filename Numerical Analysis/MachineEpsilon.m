% Finding The Machine Epsilon

clear;  % clear command window
clc;    % clear old memories to have a fresh start

% initializations ========================================================
ep = 1; % initialize the machine epsilon

% computations ===========================================================
while ep+1>1    % while ep is capable of increasing values perceptibly
    ep = ep/2;  % bit shift it down (divide by 2)
end
ep = ep*2;      % once it is too small to see, shift it back up by one bit
                % ie. the smallest amount the machine recognizes
display(ep);    % display the resulting value to the user