function [ leF, deltaA,lintact,gammath,lcontrol,Rglobal,iinit,idebond,iprop,lecrit,vcrackprop] = lefatigue( tau0, gamma0, Sinf, dSL, T, tm,  GIIth, GIIctm, G, gammaN, C1, m, lambda,  DeltaN, N, R, Lout,C,Ef,level,rf )
%Calculation of the effective recovery length increase due to the fatigue
%cycles for a specific value of remote stress
%Marco Alves
global DX
%initialize variables------------------------------------------------------
d=zeros(1,N+1);
tauD=zeros(1,N+1);
gamma=zeros(1,N+1);
lcz=zeros(1,N);
deltaA=zeros(1,N);
erate=zeros(1,N);
dAdN=zeros(1,N);
dSL=2*Sinf;
iinit=0;
idebond=0;
iprop=0;
vcrackprop=0;
%%%Calculate true stress ratio (including both effect on the fibres and matrix)
%--------------------------------------------------------------------------

%Calculate initial parameters---------------------------------------
l=xDsigma(dSL,lambda,T,tau0); % first effective recovery length
le=zeros(1,N);
lcontrol=zeros(1,N);
lintact=zeros(1,N);
le(1)=l;
lintact(1)=le(1);
tauD(1)=ftau(le(1),tau0,lambda); % Shear stress in the critical point at N=1
gammath=gammathreshold( gammaN, tau0, GIIth, GIIctm); % Calculates gamma corresponding to GIIth
Rglobal=0.1;
%1st cycle iterations
gamma(1)=fgamma(le(1),gamma0, tau0, G, lambda ); % Gamma in the critical point
 if gamma(1) < gammath % no damage can occur due to fatigue
     tauD(2:end)=tauD(1); % shear stress in the critical point
     le(2:end)=le(1);
     lintact(2:end)=le(1);
     lecrit=xtau(0,lambda,tau0);
 else
    dstatic=gamma(1)/gammaN;
    d(1)=dstatic;
    lcz1=acos((tau0^2 - 2*GIIth*abs(G))^(1/2)/tau0)/lambda;
    lecrit=xtau(0,lambda,tau0);
    %Cycle Loop------------------------------------------------------------
    for i=2:N
        d0=fdamage(gamma(i-1), (C*(le(i-1))), C1, tau0, gammaN, GIIctm, m, R )*DeltaN(i); %Damage in the critical point
        d(i)=d(i-1)+d0;
        tauD(i)=tau0*(1-d(i)); % Damaged shear stress in the critical point
        %crack propagation cylces-------------------------------------------------
        if tauD(i) < 0
%             le(i)=lecrit;
            %le(i)=xtau(tauD(i),lambda,tau0);
            le(i)=lecrit;
            lcrackprop=le(i);
            d(i)=1;
            gamma(i:end)=fgamma(lecrit,gamma0, tau0, G, lambda );
            lintact(i:end)=lecrit;
            deltaA(i)=lcrackprop-lecrit;
%           vcrackprop=(1/lambda*(acos(-fdamage(gammaN,(C*(lecrit-lcz1)),C1,tau0,gammaN,GIIctm,m,R))-pi/2));%crack propagation speed [mm/cycle]   
%           vcrackprop=fdamage(gammaN,(lecrit-lcz1),C1,tau0,gammaN,GIIctm,m,R)

%           vcrackprop=(le(i)-le(i-1))/DeltaN(i);
            vcrackprop=(le(i)-le(1))/(i*DeltaN(i));
            GIIselfSimilar = 1/(2*C)*Ef*(Sinf/Ef)^2*2^(level-1)*pi*rf^2;
            Gratio=GIIselfSimilar/(GIIctm*tm);
            vcrackprop=C1/C*(Gratio)^m;
            
%             if vcrackprop < vcrackbut1
%                 vcrackprop = 2*vcrackbut1-vcrackbut2;
%             end
            for k=i+1:N
                tauD(k)=0;
                le(k)=le(k-1)+vcrackprop*DeltaN(k);
                deltaA(k)=deltaA(k-1)+vcrackprop*DeltaN(k);
                d(k)=1;
                iprop=iprop+DeltaN(k);
            end
            break 
        end
     le(i)=xtau(tauD(i),lambda,tau0);
     lintact(i)=le(i);
     gamma(i)=fgamma(le(i),gamma0, tau0, G, lambda );
     iinit=iinit+DeltaN(i);
    end
 end
 
 %Correct the values of deltaA and lintact since the sum of both can never
 %be bigger than Lout
boundary=(Lout/2)-(xtau(0,lambda,tau0));
 for i=1:N
     if deltaA(i) > boundary 
         le(i:end)=Lout/2;
         lintact(i)=lintact(i)-(deltaA(i)-boundary);
         if lintact(i) < 0
             lintact(i:end)=0;
             deltaA(i:end)=Lout/2;
             for o=i:N
                 idebond=idebond+DeltaN(o);
             end
             break
         end
     end 
 end
 
%Multiply by 2 because of the symmetry considered in the analysis
leF=2*le;
deltaA=2*deltaA;
lintact=2*lintact;

%Definition of the control length
lcontrol=2*leF;

for i=1:N
    if lcontrol(i)>Lout
        lcontrol(i)=Lout;
    end
end

%Energy release rate calcualtion
% for i=1:N
%     erate(i)=ERR(gamma(i),gammaN,tau0);
%     dAdN(i)=fdAdN( gamma(i), C1, tau0, gammaN, GIIctm, m, R ); 
% end


end

