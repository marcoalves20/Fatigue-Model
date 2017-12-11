function [ critstress ] = criticalstress( tau0, lambda, Tf  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

critstress=tau0./(lambda.*Tf);

end

