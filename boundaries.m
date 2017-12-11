function [] = boundaries(n,A,C,Ccorner,Cedge)
%runs bundle strength model for free-boundaries effects
%called by bundleX function (T=41, freebounds=1)

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global DX 
global Lref 
%Internal
global nX 
%Output
global X Lout
%Input from NORMAL SLB
global lnSur lnSlr lnSkr
global lnSue lnSke le
%Output for SUPER-BUNDLE SLB
global lnSurS lnSlrS lnSkrS
global FuoS XavgoS CoVoS lnSuoS

%%%Runs model with CORNER SLB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preallocating matrices
lnSurC=zeros(nX,n-1);
lnSlrC=zeros(nX,n-1);
lnSkrC=zeros(nX,n-1);
lnSueC=zeros(nX,n-1);
lnSkeC=zeros(nX,n-1);
leC=zeros(nX,n-1);

%Level zero
[lnSurC(:,1),lnSlrC(:,1),lnSkrC(:,1),leC(:,1),lnSueC(:,1),lnSkeC(:,1)]=...
    levelzero(A(1),Ccorner(1));

%Level 1:n-2
for i=1:n-2
    j=i+1;
    [lnSurC(:,j),lnSlrC(:,j),lnSkrC(:,j),leC(:,j),lnSueC(:,j),lnSkeC(:,j)]=...
        leveli(lnSue(:,j-1),lnSke(:,j-1),le(:,j-1),A(j),Ccorner(j));
end

%%%Runs model with EDGE SLB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preallocating matrices
lnSurE=zeros(nX,n-1);
lnSlrE=zeros(nX,n-1);
lnSkrE=zeros(nX,n-1);
lnSueE=zeros(nX,n-1);
lnSkeE=zeros(nX,n-1);
leE=zeros(nX,n-1);

%Level zero
[lnSurE(:,1),lnSlrE(:,1),lnSkrE(:,1),leE(:,1),lnSueE(:,1),lnSkeE(:,1)]=...
    levelzero(A(1),Cedge(1));

%Level 1:n-1
for i=1:n-1
    j=i+1;
    [lnSurE(:,j),lnSlrE(:,j),lnSkrE(:,j),leE(:,j),lnSueE(:,j),lnSkeE(:,j)]=...
        leveli(lnSueC(:,j-1),lnSkeC(:,j-1),leC(:,j-1),A(j),Cedge(j));
end


%%%Runs model with SUPER-BUNDLE SLB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preallocating matrices
lnSurS=lnSur;
lnSlrS=lnSlr;
lnSkrS=lnSkr;

%Level 1:n
for i=1:n
    j=i+1;
    [lnSurS(:,j),lnSlrS(:,j),lnSkrS(:,j)]=...
        leveli(lnSueE(:,j-1),lnSkeE(:,j-1),leE(:,j-1),A(j),C(j));
end

%%%Post-process results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bundle strength distribution for output length
lnSuoS=Lout/Lref*lnSurS;
FuoS=1-exp(lnSuoS);

%Average bundle strength for output length
intFuoDXS=trapz(FuoS)*DX;
XavgoS=X(end)-intFuoDXS;

%CoV of bundle strength for output length
intXFuoDXS=zeros(1,n+1);
for i=1:n+1
    intXFuoDXS(1,i)=trapz(X.*FuoS(:,i))*DX;
end

CoVoS=sqrt((X(end)-XavgoS).^2-2.*intXFuoDXS+2.*XavgoS.*intFuoDXS)./XavgoS;


end




