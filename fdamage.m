function [ dDdN ] = fdamage( Gamma, Acz, C1, tau0, gammaN, GIIctm, m, R )

% dDdN=(1/Acz)*C1*(((tau0/2)*(gammaN-(gammaN-Gamma)^2/gammaN)*(1-R^2))/GIIctm)^m;

deltaG=tau0/2*(gammaN-((gammaN-Gamma)^2/gammaN))*(1-R^2);
dDdN=(1/Acz)*C1*(deltaG/GIIctm)^m;

end

