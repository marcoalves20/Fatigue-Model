function [gamma] = fgamma( x,gamma0, tau0, G, lambda )

gamma = gamma0+(tau0/abs(G))*(1-cos(lambda*x));
