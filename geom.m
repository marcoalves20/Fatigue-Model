function[A,C,nf,Ccorner,Cedge,Tf,tm]=geom(n,c,T,Df,Vf)
%calculates geometry of all bundle levels
%called by bundleX

%global AMu

%Calculates level-0 geometry
Af=pi/4*Df^2;
Cf=pi*Df;

%Calculates geometric parameters for each level
nf=c.^(0:1:n);
A=nf.*Af;

%Calculates arrangement geometry
if T==41||T==43||T==45 % quadrangular case
    lQ=sqrt(pi)/(2*sqrt(Vf))*Df;
    sQ=(sqrt(pi)/(2*sqrt(Vf))-1)*Df;
    RnfQ=sqrt(nf);
elseif T==61 % hexagonal case
    lH=sqrt(pi/(6*sqrt(3)*Vf))*Df;
    sH=(sqrt(pi/(2*sqrt(3)*Vf))-1)*Df;
    RnfH=nf.^(log(3)/log(7));
end

%Calculates circumference
if T==43 % quadrangular, matrix
    C=4*RnfQ.*lQ;
    Ccorner=0;
    Cedge=0;
%    AMu=1;
elseif T==41 % quadrangular, interfacial
    C=3*Cf+4*((RnfQ-1)*sQ+(RnfQ-2)*Cf/2);
    Ccorner=Cf/4 + 2*(RnfQ-1/2)*(sQ+Cf/2);
    Cedge= sqrt(2)*RnfQ*(sQ+Cf/2);
%    AMu=(sQ+Df)/(sQ+Cf/2);
elseif T==45 % quadrangular, shortest
    C=Cf+4*(RnfQ-1)*lQ;
    Ccorner=0;
    Cedge=0;
%    AMu=1;
elseif T==61 % hexagonal, interfacial
    C=3*(RnfH-1)*sH+(3*RnfH-1)*Cf/2; 
    Ccorner=0;
    Cedge=0;
%    AMu=(sH+Df)/(sH+Cf/2);
end

Tf=A./C;
tm=sQ;
% Tf1=Tf;

%Shift the values of Tf to the right Tf=[0,0,1,2,...i]
% Tf1=Tf;
% Tf(1:2)=Tf1(1);
% Tf(3:length(A)+1)=Tf1(2:end);



end







