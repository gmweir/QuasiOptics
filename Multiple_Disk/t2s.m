function [ s ] = t2s( t )
%T2S Summary of this function goes here
%    By Fabian Wilde (IPP has no RF Toolbox)
    s = zeros(2,2);
    s(1,1) = t(1,2)/t(2,2);
    s(1,2) = t(1,1) - ((t(1,2)*t(2,1))/(t(2,2)));
    s(2,1) = 1/t(2,2);
    s(2,2) = -t(2,1)/t(2,2);
end