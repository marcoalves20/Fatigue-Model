function[] = microcomp(kXavg,kCoV)
%runs bundle strength model for validation at micro-scale
%Data from Beyerlein and Phoenix (1996a) and Kazanci (2004) 

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global Xmax DX
global Lin Xref m
global Lref Lout
global Tsl T freebounds
global Df Vf
global k

%Outputs
global lnSuoS X
global lnlnSuoSout


%%%Fixed inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xmax=50000; DX=1;
Lin=10; Lref=1; Lout=10; 
T=41; freebounds=1;
k=2;

%%%Variable inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fibre
XavgS=4116; CoVS=0.238;
XavgI=4872; CoVI=0.213;

%Resin
TslS=46.6; 
TslF=10.3;

%Bundle
DfQ=0.00685; VfQ=0.70; nfbQ=4;
DfH=0.00563; VfH=0.56; nfbH=7;


%%%Preliminary calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fibre strength distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Scaling the MLM fibre properties
XavgS=kXavg*XavgS; CoVS=kCoV*CoVS;
XavgI=kXavg*XavgI; CoVI=kCoV*CoVI;

%Finding distribution parameters
%S fibres
mS=fzero(@(mS) sqrt(gamma(1+2/mS)/(gamma(1+1/mS)^2)-1)-CoVS, 1.2/CoVS);
XinS=XavgS/(gamma(1+1/mS));
XrefS=XinS*((Lref/Lin)^(-1/mS)); 
%I fibres
mI=fzero(@(mI) sqrt(gamma(1+2/mI)/(gamma(1+1/mI)^2)-1)-CoVI, 1.2/CoVI);
XinI=XavgI/(gamma(1+1/mI));
XrefI=XinI*((Lref/Lin)^(-1/mI)); 

%%%Defining vectors with variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ProVec=[SSQ,   SFQ,    ISH,    IFH]    
XrefVec=[XrefS,	XrefS,  XrefI,	XrefI];
mVec=   [mS,	mS,     mI,     mI];
TslVec=	[TslS,	TslF,	TslS,	TslF];
DfVec=	[DfQ,	DfQ,	DfH,	DfH];
VfVec=	[VfQ,	VfQ,	VfH,	VfH];
nfbVec=	[nfbQ,	nfbQ,	nfbH,	nfbH];


%Preallocating matrices
lnlnSuoSnf=zeros(2000,4);
lnlnSuoSi0=zeros(2000,2);



%%%Running the model for each type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for itype=1:4
    
    %Clears variables
    clearvars -global -except ...
        Xmax DX Lin Lref Lout T freebounds k ...
        XrefVec mVec TslVec DfVec VfVec nfbVec...
        Xref m Tsl Df Vf nfb ...
        itype ...
        X lnSuoS lnlnSuoSnf lnlnSuoSi0 lnlnSuoSout
    
    % Defines properties for this specific type
    Xref=XrefVec(itype); m=mVec(itype); Tsl=TslVec(itype);
    Df=DfVec(itype); Vf=VfVec(itype); nfb=nfbVec(itype);
    
    % Runs bundle model
    bundleX(4);
    lnlnSuoSnf(:,itype)=interp1([1:1:4]',log(abs(lnSuoS(1:10:20000,2:1:5)))',log2(nfb),'spline')';
    if itype==1
        lnlnSuoSi0(:,1)=log(abs(lnSuoS(1:10:20000,1)));
    elseif itype==3
            lnlnSuoSi0(:,2)=log(abs(lnSuoS(1:10:20000,1)));
    end        
end

lnlnSuoSout=[log(X(1:10:20000)),lnlnSuoSi0,lnlnSuoSnf];

end

