function [ x ] = xDsigma(knownsigma,lambda,T,tau0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x=1/lambda*asin(lambda*T*knownsigma/(2*tau0));
end

