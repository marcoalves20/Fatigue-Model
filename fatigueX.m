function [] = fatigueX(n,varinow)
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
global N DeltaN
global sigmath

%%%Input variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fatigue Cycles Defenition
%Variable jump 0cycle 

DeltaN=zeros(5000,1);
DeltaN(1:5000)=10;
% DeltaN(10001:15000)=500;
% DeltaN(2001:3000)=100;
% DeltaN(3001:3500)=1000; 
% DeltaN(3501:4000)=10000; 
%--------------------

N=length(DeltaN); % Number of iterations
R=0.1; % Fatigue Stress Ratio

%Numerical variables
Xmax=40000; DX=100;

%Input fibre strength distribution
% Lin=5; Xavg=8500; CoV=0.25; %Original Parametric study

Lin=20; Xavg=6200; CoV=0.25; %Original

%Interfacial shear strengthnn
Tsl=70;

%Geometry and composite
Df=0.0071; Vf=0.6; % T700 CF 

%Shear lag boundary:
%4-quadrangular, 6-hexagonal;
%1-interface, 3-matrix, 5-shortest.
T=41; freebounds=0;

%Lengths
%T700G Lref = 66
%T700S Lref = 50
%T700 Yuaxin Zhou Lref = 8
%Reference length (unitary)
Lref=1; 
% Output length
Lout=180; 

%Stress concentrations
k=varinow;

%%%Preliminary calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fibre strength distribution
% T700G==========
% m=8.44;       %=
% Xref=5054;    %=
%===============
% T700S=========
% m=4.8;       %=
% Xref==========
% m=4.8; Original Parametric Study  
m=5;

% m=fzero(@(m) sqrt(gamma(1+2/m)/(gamma(1+1/m)^2)-1)-CoV, 1.2/CoV)
Xin=Xavg/(gamma(1+1/m));
Xref=Xin*((Lref/Lin)^(-1/m)); 

%%%Calculates strength distributions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bundleX(n,DeltaN,N,R,Df);
toc
end