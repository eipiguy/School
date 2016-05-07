function [J,gx]= g(x)
gx=[sqrt(5-x(2).^(2));x(1).^(2)-1];
J=[0,(-x(2)./(sqrt(5-x(2).^(2))));2*x(1),0];
end

