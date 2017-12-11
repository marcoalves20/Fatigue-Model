function [ leF, deltaA,lintact,gammath,lcontrol,Rglobal,iinit,idebond,iprop,lecrit,vcrackprop] = lefatigue( tau0, gamma0, Sinf, dSL, T, tm,  GIIth, GIIctm, G, gammaN, C1, m, lambda,  DeltaN, N, R, Lout,C,Ef,level,rf,A )
%Calculation of the effective recovery length increase due to the fatigue
%cycles for a specific value of remote stress
%Marco Alves
global DX
%initialize variables------------------------------------------------------
res=10; % resolution of domain discritization
d=zeros(1,N+1);
tauD=zeros(res,N+1);
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
lecrit=0;
%%%Calculate true stress ratio (including both effect on the fibres and matrix)
%--------------------------------------------------------------------------

%Calculate initial parameters---------------------------------------
l=xDsigma(dSL,lambda,T,tau0); % first effective recovery length
le=zeros(1,N);
lcontrol=zeros(1,N);
lintact=zeros(1,N);
le(1)=l;
x = linspace(0,le(1),res);%Discritize domain 
resolution = x(2)-x(1);
lintact(1)=le(1);
tauD(:,1)=ftau(x,tau0,lambda); % Shear stress in the critical point at N=1
damage(:,1)=1-tauD(:,1)/tau0;
area(1)=resolution*trapz(tauD(:,1))*C/A;
gammath=gammathreshold( gammaN, tau0, GIIth, GIIctm); %Calculates gamma corresponding to GIIth
Rglobal=0.1;

%1st cycle iterations
gamma(1)=fgamma(le(1),gamma0, tau0, G, lambda ); % Gamma in the critical point

 if gamma(1) < gammath % no damage can occur due to fatigue
     le(2:end)=le(1);
 else
    %Cycle Loop------------------------------------------------------------
    for i=2:N
        for j=1:length(x)
            damage(j,i) = fdamagetau(tauD(j,i-1), (C*le(i-1)), C1, tau0, gammaN, GIIctm, m, R )*DeltaN(i);
        end
        aux = sum(damage,2);
        tauD(:,i) = (1-aux).*tau0;
        %crack propagation cylces------------------------------------------
        if tauD(length(x),i)<0
            area(i) = resolution*trapz(tauD(:,i))*C/A;
            le(i) = le(1)*(area(1)/area(i));
            %Interpolation to get the exact value of le for tau=0(@crittical point)
            interpY=[le(i-1) le(i)];
            interpX=[tauD(length(x),i-1) tauD(length(x),i)];
            le(i) = interp1(interpX,interpY,0);
            lecrit=le(i);
            lintact(i:end)=lecrit;
            deltaA(i)=0;
            % Crack propagation definition
            vcrackprop=(le(i)-le(1))/(i*DeltaN(i));
            GIIselfSimilar = 1/(2*C)*Ef*(Sinf/Ef)^2*2^(level-1)*pi*rf^2;
            Gratio=GIIselfSimilar/(GIIctm*tm);
            vcrackprop=C1/C*(Gratio)^m;
            for k=i+1:N
                le(k)=le(k-1)+vcrackprop*DeltaN(k);
                deltaA(k)=deltaA(k-1)+vcrackprop*DeltaN(k);
            end
            break 
        end
     area(i) = resolution*trapz(tauD(:,i))*C/A;
     le(i) = le(1)*(area(1)/area(i));
     lintact(i)=le(i);
    end
 end
 
 %Correct the values of deltaA and lintact since the sum of both can never
 %be bigger than Lout
boundary=(Lout/2)-(lecrit);
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

