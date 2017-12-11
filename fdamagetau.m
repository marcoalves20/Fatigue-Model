function [ dDdN ] = fdamage( tau, Acz, C1, tau0, gammaN, GIIctm, m, R )

%FDDDN
gammaN = 2*GIIctm/tau0;
Gdm = tau0/gammaN;
gammai = (tau0-tau)/Gdm;
% 
% dDdN=(1/Lcz)*C1*(((tau0/2)*(gammaN-(gammaN-tau)^2/gammaN))/GIIctm*(1-R^2))^m;
deltaG = tau0/2 * (gammaN - ((gammaN-gammai)^2)/gammaN) * (1-R^2);
dDdN=(1/Acz)*C1*(deltaG/GIIctm)^m;

end

