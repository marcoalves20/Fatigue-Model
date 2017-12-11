function [lnSurI,lnSucI] =...
    leveli(lnSueO,lnSke1O,lnSke2O,leO,lcontrolO,deltaAO,lintactO)
%xxI -> related to level i+1
%xxO -> related to level i

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%Internal
global nX N
%Output

%==============================Scalling Law================================
%auxiliar parameters - all lengths are related to level-[i]
a=(lcontrolO./leO).*lnSueO;
b=(deltaAO./leO).*lnSke2O;
c=(lintactO./leO).*lnSke1O;
d=((lcontrolO-leO)./leO).*lnSueO;

%--------------------------------------------------------------------------
%Suvival probability of level-[i+1]
lnSucI=2*a+log(1-2*exp(c+d+b-a)+2*exp(c+d+b-2*a));
% lnSucI=2*a+c+d+log(exp(-c-d)-2*exp(b-a)+2*exp(b-2*a));
% lnSucI(lnSucI==Inf)=2*a(lnSucI==Inf)+c(lnSucI==Inf)+d(lnSucI==Inf)+log(exp(-c(lnSucI==Inf)-d(lnSucI==Inf))-2*exp(b(lnSucI==Inf)-a(lnSucI==Inf))+2*exp(b(lnSucI==Inf)-2*a(lnSucI==Inf)));
% lnSucI(lnSucI==Inf)=2*a(lnSucI==Inf)+c(lnSucI==Inf)+log(exp(-c(lnSucI==Inf))+2*exp(b(lnSucI==Inf)+d(lnSucI==Inf)-2*a(lnSucI==Inf))-2*exp(b(lnSucI==Inf)+d(lnSucI==Inf)-a(lnSucI==Inf)));
% lnSucI(lnSucI==Inf)=2*a(lnSucI==Inf)+log(1-2*exp(c(lnSucI==Inf)+d(lnSucI==Inf)+b(lnSucI==Inf)-a(lnSucI==Inf))+2*exp(c(lnSucI==Inf)+d(lnSucI==Inf)+b(lnSucI==Inf)-2*a(lnSucI==Inf)));
lnSucI(1)=0;

%Scalling back to the reference length-------------------------------------

% lnSurI=zeros(nX,1); % I changed this part here.. Change back to nX if
% looping in cycles and not in sigma.
lnSurI=zeros(nX,1);
lnSurI(:,1)=1./(lcontrolO).*lnSucI;
lnSurI(1)=0;
lnSurI(lnSurI>0)=-lnSurI(lnSurI>0);
%==========================================================================
lnSurI = sort(lnSurI,'descend');
end

