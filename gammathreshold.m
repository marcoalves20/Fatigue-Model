function [ gammath ] = gammathreshold( gammaN, tau0,GIIth,GIIc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%   gammath=gammaN-sqrt(gammaN*(gammaN-2*GIIth/tau0));
% gammath=-((-tau0*gammaN+sqrt(gammaN^2*tau0^2-2*GIIth*gammaN*tau0))/tau0);
gammath=-(sqrt(2)*sqrt((GIIc-GIIth)*tau0*gammaN)-tau0*gammaN)/tau0;
end

