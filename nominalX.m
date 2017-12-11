function [] = nominalX(n,varinow)
%defines inputs for strength model, and runs bundleX
tic
%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global Xmax DX
global Lin Xavg CoV m Xin
global Lref Xref Lout
global Tsl T freebounds
global Df Vf
global k
global N

%%%Input variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fatigue Cycles Defenition
DeltaN=zeros(10,1);
DeltaN(1:10)=1;
%--------------------
N=length(DeltaN); % Number of iterations % Number of iterations. Total N of cycles = N*DeltaN; % Number os cycles per iteration during stable crack propagation
R=0; % Fatigue Stress Ratio

%Numerical variables
Xmax=50000; DX=10;

%Input fibre strength distribution
Lin=10; Xavg=4500; CoV=0.25; 

%Interfacial shear strength
Tsl=70;

%Geometry and composite
Df=0.005; Vf=0.6;

%Shear lag boundary:
%4-quadrangular, 6-hexagonal;
%1-interface, 3-matrix, 5-shortest.
T=41; freebounds=0;

%Lengths
Lref=1; Lout=10; 

%Stress concentrations
k=varinow;

%%%Preliminary calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fibre strength distribution
m=fzero(@(m) sqrt(gamma(1+2/m)/(gamma(1+1/m)^2)-1)-CoV, 1.2/CoV);
Xin=Xavg/(gamma(1+1/m));
Xref=Xin*((Lref/Lin)^(-1/m)); 

%%%Calculates strength distributions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bundleX(n,DeltaN,N,R);
toc
end