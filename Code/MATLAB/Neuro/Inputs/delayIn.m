function [ InV ] = delayIn( t, V, d )
    m = length(V);
    InV = zeros(m,1);
    for i=1:m
        if t < d(i)
            InV(i) = 0;
        else
            InV(i) = V(i);
        end
    end
end