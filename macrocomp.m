function [] = macrocomp(n)
%runs bundle strength model for validation at macro-scale
%Data from Okabe and Takeda (2002)

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global Xmax DX
global Lin m Xin
global Lref Xref Lout
global Tsl T freebounds
global Df Vf
global k
%Outputs (confidence bounds)
global X alpha lnSuo boundso

%%%Input variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Numerical variables
Xmax=5000; DX=100;

%Input fibre strength distribution
Lin=50; m=3.8; Xin=3570;

%Interfacial shear strength
Tsl=52.4*1.0;

%Geometry and composite
Df=0.005; Vf=0.6;

%Shear lag boundary:
%4-quadrangular, 6-hexagonal;
%1-interface, 3-matrix, 5-shortest.
T=41; freebounds=0;

%Lengths
Lref=1; Lout=10; 

%Stress concentrations
k=2;

%%%Preliminary calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xref=Xin*((Lref/Lin)^(-1/m)); 

bundleX(n)

%Boundaries for confidence intervals at alpha,1-alpha

alpha=(0.05:0.05:0.95)';
for i=0:n
    boundso(:,i+1)=interp1(lnSuo(:,i+1),X,log(1-alpha));
end

end