function [ Rglobal ] = trueratio( R,Sinf,T,gamma0, tau0, G, lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global nX DX

Smax=Sinf;
Smin=Smax*R;

lmax=xDsigma(2*Smax,lambda,T,tau0);
lmin=xDsigma(2*Smin,lambda,T,tau0);

GammaMax=fgamma(lmax,gamma0, tau0, G, lambda);
GammaMin=fgamma(lmin,gamma0, tau0, G, lambda);

Rglobal=GammaMin/GammaMax;
    

end

