function [ ER ] = ERR( gamma, gammaN, tau0 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
ER=(tau0/2)*(gammaN-((gammaN-gamma)^2/gammaN));

end

