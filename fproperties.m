function [ tau0, gamma0, dSL, GIIth, GIIctm, G, gammaN, C1, mp,  lambda,E ] = fproperties( X, Tf,tm,Tsl )
%Marco Alves

%==============Fiber Young Modulus==============
E=240e3;%[MPa] % T700 CF
%===============================================

%==============B.C @ xi=0==============
dS0=0; %[Pa];
tau0=Tsl;%[MPa];
gamma0=0;
%============BC @ xi/le/2=============
dSL=2*X;%[MPa];
%======================================

%==========Matrix Properties===========
GIIc=1;%[kJ/m^2];
GIIctm=GIIc/tm;
gammaN=(GIIctm)*2/(tau0);
% GIIthresh=0.02;
GIIthresh=0.02;
GIIth=GIIthresh/tm;
G=-(tau0/gammaN);%[MPa];
%======================================

%==========Fatigue Parameters==========
C1=4.5e-5;
% C1=4.5e-3;

mp=1.5;
% mp=3.5;
x0=0;
lambda=sqrt(2*abs(G)./(Tf*tm*E));
%======================================
end

