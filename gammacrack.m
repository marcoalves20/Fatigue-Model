function [ gamma ] = gammacrack( x, x0, gammaN, sinf, tm, E )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

gamma=gammaN+sinf/(tm*E)*(x-x0);
end

