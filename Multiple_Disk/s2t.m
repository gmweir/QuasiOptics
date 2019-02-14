function [ t ] = s2t( s )
%S2T Summary of this function goes here
%   By Fabian Wilde (IPP has no RF Toolbox)
    t = zeros(2,2);
    t(1,1) = s(1,2) - ((s(1,1)*s(2,2))/(s(2,1)));
    t(1,2) = s(1,1)/s(2,1);
    t(2,1) = -s(2,2)/s(2,1);
    t(2,2) = 1/s(2,1);
end