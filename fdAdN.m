function [ dDdN ] = fdAdN( Gamma, C1, tau0, gammaN, GIIctm, m, R )

%FDDDN

dDdN=C1*(((tau0/2)*(gammaN-(gammaN-Gamma)^2/gammaN))/GIIctm*(1-R^2))^m;

end

