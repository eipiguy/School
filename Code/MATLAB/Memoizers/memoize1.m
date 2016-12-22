function f = memoize1(F)
%% memoize1: adds a cache to a one-arg function F, inputs testable with ==
%   [ f ] = memoize1( F )
%       Returns a cached/memoized version of the input function,
%       where inputs and output are stored in arrays,
%       and equivalences for each new input are checked with ==
%##########################################################################
%% Input:
%   F = function to cache
%##########################################################################
%% Output:
%   f = cached version of F
%##########################################################################

% Memory cache initialization
inMem = [];
outMem = [];

% Start the cached function as the input
f = @inner;
    % Then rewrite the inner function to store into memory, 
    % and check memory before returning:
    function out = inner(in)
        out = zeros(size(in));  % preallocate output
        
        % Find indices of inputs that are already in memory.
        [tf,loc] = ismember(in,inMem);         
        ft = ~tf;   % Then the indices of new inputs the other ones.
        
        out(ft) = F(in(ft));    % Caclculate outputs for new inputs.
        
        % Place new values in memory.
        inMem = [inMem in(ft(:).')];
        outMem = [outMem reshape(out(ft),1,[])];
        
        out(tf) = y(loc(tf));  % Fill in the memorized output values.
    end
end