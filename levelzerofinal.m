function [lnSur,lnSue,lnSke1,lnSke2] =...
    levelzero(AI,CI,leF,lcontrol,deltaA,lintact)
%runs bundle strength model for level zero (single fibre, Weibull)
%called by master file bundleX

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global Tsl
global Lref
%Internal
global m Xref
global Ck nK nX DX
%Output
global X k

%%%For the reference length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates survival functions
lnSkr1=zeros(nX,1);
lnSkr2=zeros(nX,1);

% lnSur=-((X/Xref).^m);
%=========================================
lnSlr=-1/(m+1)*((X/Xref).^m);

lnSur=-((X/Xref).^m);
lnSur(1)=0;
%=========================================
lnSlrI=cumtrapz(lnSur(1:nX,1))*DX./X;
        lnSlrI(1)=0;
        if mod(k,1)==0
            %If k is integer, takes every other k-value:
            lnSkr1(1:nK,1)=(k*lnSlrI(1:k:end)-lnSlrI(1:nK))./(k-1);
        else
            %If k is not integer (interpolation required):
            %Interpolates lnSlr for k*X (non-integer k):
            lnSlrKI(1:nK,1)=interp1(X,lnSlrI,XK,'linear');
            %uses lnSlr interpolated @ k*X (lnSlrK)
            lnSkrI1(1:nK,1)=(k*lnSlrKI-lnSlrI(1:nK))./(k-1);
        end
lnSkr1(nK+1:nX,1)=Ck*lnSur(nK+1:nX,1);
%=========================================
% lnSkr1=-Ck*((X/Xref).^m);


lnSkr2(1:nK,1)=lnSur(1:k:end);
lnSkr2(nK+1:end,1)=k^m*lnSur(nK+1:end);
lnSkr1(1)=0;
lnSkr2(1)=0;

%%%For effective recovery length:

%Calculates survival functions
% lnSue=leF.*lnSur;
% lnSke1=lintact.*lnSkr1;
% lnSke2=deltaA.*lnSkr2;

lnSue=leF.*lnSur;
lnSke1=leF.*lnSkr1;
lnSke2=leF.*lnSkr2;

 end

