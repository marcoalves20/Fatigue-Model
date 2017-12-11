function [lnSurI,lnSlrI,lnSkrI,lnSueI,lnSkeI1,lnSkeI2] =...
    leveli(lnSueO,lnSke1O,lnSke2O,AI,CI,leO,leI,lcontrolO,lcontrolI,deltaAO,deltaAI,lintactO,lintactI,nXcritical,lambda,Tf,tau0)
%xxI -> related to level i+1
%xxO -> related to level i

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global DX
global Tsl
global Lref Lout
global k m
%Internal
global Ck 
global nX nK XK 
%Output
global X Xref 

%==============================Scalling Law================================
for j=1:nX
    sigmaK(j)=k*X(j);
end

a=(lcontrolO./leO).*lnSueO;
b=(deltaAO./leO).*lnSke2O;
c=(lintactO./leO).*lnSke1O;
d=((lcontrolO-leO)./leO).*lnSueO;

lnSucI=2*a+c+d+log(exp(-c-d)-2*exp(b-a)+2*exp(b-2*a));
lnSucI(lnSucI==Inf)=2*a(lnSucI==Inf)+c(lnSucI==Inf)+log(exp(-c(lnSucI==Inf))+2*exp(b(lnSucI==Inf)+d(lnSucI==Inf)-2*a(lnSucI==Inf))-2*exp(b(lnSucI==Inf)+d(lnSucI==Inf)-a(lnSucI==Inf)));
lnSucI(lnSucI==Inf)=2*a(lnSucI==Inf)+log(1-2*exp(c(lnSucI==Inf)+d(lnSucI==Inf)+b(lnSucI==Inf)-a(lnSucI==Inf))+2*exp(c(lnSucI==Inf)+d(lnSucI==Inf)+b(lnSucI==Inf)-2*a(lnSucI==Inf)));

lnSucI(1)=0;
%Scalling back to the reference length-------------------------------------
lnSurI=zeros(nX,1);
lnSurI(:,1)=1./(lcontrolO).*lnSucI;
lnSurI(1)=0;
%-------------------------------------------------------------------------

lnSkrI1=zeros(nX,1);
lnSkrI2=zeros(nX,1);
lnSkeI1=zeros(nX,1);
lnSkeI2=zeros(nX,1);
%Calculates SkrI1 up to X and then it is limited to nXcrit in SkrI1
lnSkrI1=zeros(nX,1);
 lnSlrI=cumtrapz(lnSurI(1:nX))*DX./X;
 lnSlrI(1)=0;
    if mod(k,1)==0
        %If k is integer, takes every other k-value:
        for i = 1:nK-1
             lnSkrI1(i,1)=2*Tf/(tau0*(leI(i)/2))*trapz(lnSurI(i:(2*i))/sqrt(4-(lambda^2*Tf^2*(k*X(i)-2*sigmaK(i))^2/tau0^2)))*DX; 
        end
    else
        %If k is not integer (interpolation required):   
        %Interpolates lnSlr for k*X (non-integer k):
        lnSlrKI(1:nK,1)=interp1(X,lnSlrI,XK,'linear');
        %uses lnSlr interpolated @ k*X (lnSlrK)
        lnSkrI1(1:nK,1)=(k*lnSlrKI-lnSlrI(1:nK))./(k-1);
    end
    lnSkrI1(nK:nX,1)=Ck*lnSurI(nK:nX);
    lnSkrI1(1)=0;

    %----------------------------------------------------------------------  
    %Survival probability calculated for the stress plateau
    %Modification 21/03/17-------------------------------------------------
    lnSkrI2(1:nK)=lnSurI(1:k:end);
    lnSkrI2(nK+1:end)=k^m*lnSurI(nK+1:end);
    lnSkrI2(1)=0;
    %----------------------------------------------------------------------
    lnSkrI=lnSkrI1;
    %%%For effective recovery length%%%%%%%%%%%%%%%%%%
    %Calculates survival functions    
    lnSkeI1=lintactI.*lnSkrI1;
    lnSkeI2=deltaAI.*lnSkrI2;
    lnSueI=(leI).*lnSurI;

end

