function [J,fx]= f(x)
fx=[sqrt(x(2)+1);sqrt(5-x(1).^(2))];
J=[0,(1./(2*sqrt(x(2)+1)));(-x(1)./sqrt(5-x(1).^(2))),0];
end

