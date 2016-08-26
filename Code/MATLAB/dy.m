function [ dydt ] = dy( t,y0,varargin )
    dydt = -0.5*(eye(2)*y0');
end

