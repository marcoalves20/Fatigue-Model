function [lnSurI,lnSlrI,lnSkrI,lnSueI,lnSkeI,lnSue1] =...
    levelzero(AI,CI,leF,lcontrol)
%runs bundle strength model for level zero (single fibre, Weibull)
%called by master file bundleX

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global Tsl
global Lref
%Internal
global m Xref
global Ck
%Output
global X 

%%%For the reference length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates survival functions
    lnSurI=-((X/Xref).^m);
    lnSlrI=-1/(m+1)*((X/Xref).^m);
    lnSkrI=-Ck*((X/Xref).^m);

%%%For effective recovery length:

%Calculates survival functions
lnSueI=leF.*lnSurI;
lnSkeI=leF.*lnSkrI;
lnSue1=(lcontrol-leF).*lnSurI;

 end

