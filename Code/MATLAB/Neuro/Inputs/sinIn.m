function [ InV ] = sinIn( t, a, f, p )
    m = length(a);
    InV = zeros(m,1);
    for i=1:m
        InV(i) = a(i)*sin( (f(i)*t) + p(i) );
    end
end