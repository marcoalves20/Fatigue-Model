function [] = bundleX(n,DeltaN,N,R,Df)
%runs bundle strength model until level n

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global Xmax DX
global m Tsl
global Lref Lout
global T freebounds
global Df Vf
global k
%Internal
global Ck
global nX nK XK
%Output
global X
%Output for NORMAL SLB
global lnSur lnSuc lnSlr lnSkr
global lnSue lnSke le leF deltaA lintact
global Fuo Xavgo CoVo lnSuo
global C A
global ER dAdN GIIctm GIIth
global sigmath
global nXcrit
global lnSke1 lnSke2
global iinit idebond iprop lecrit vcrackprop



%=========================Preliminary calculations=========================
%Stress concentrations factor
if k==1
    Ck=1;
else
    Ck=(k.^(m+1)-1)/((k-1)*(m+1));
end

%==============================Loading vector==============================
%Number of load increments
nX=ceil(Xmax/DX)+1;

%Number of load increments for Sk
nK=floor((nX-1)/k+1);

%Generates the load vector
X=transpose(linspace(0,DX*(nX-1),nX));

%Load vector for stress concentrations
XK=k*X(1:nK);

%========================Geometry of fibre bundles=========================
[A,C,nf,Ccorner,Cedge,Tf,tm]=geom(n,2,T,Df,Vf);

%=====Definition of the Cohesive law parameters and Paris law constants====
[ tau0, gamma0, dSL, GIIth, GIIctm, G, gammaN, C1, mp,  lambda,Ef ] = fproperties( X,Tf,tm,Tsl );

%==========Definition of the critical stress for each bundle level=========
critstress = criticalstress(tau0, lambda, Tf);
nXcrit=floor(critstress./DX)+1;

%============================Fatigue Cycles================================
%Preallocating matrices----------------------------------------------------
leF=zeros(nX,n+1,N);
lcontrol=zeros(nX,n+1,N);
deltaA=zeros(nX,n+1,N);
lintact=zeros(nX,n+1,N);
ER=zeros(nX,n+1,N);
dAdN=zeros(nX,n+1,N);
Rglobal=zeros(nX,n+1);
iinit=zeros(nX,n+1);
idebond=zeros(nX,n+1);
iprop=zeros(nX,n+1);
lecrit=zeros(nX,n+1);
thindex=zeros(n+1);
vcrackprop=zeros(nX,n+1);
%Calculation of the all the effective recovery lengths---------------------
for i=1:nX
    for j=1:n+1
        if i <= nXcrit(j)
            [ leF(i,j,:), deltaA(i,j,:), lintact(i,j,:),gammath,lcontrol(i,j,:),Rglobal(i,j),iinit(i,j),idebond(i,j),iprop(i,j),lecrit(i,j), vcrackprop(i,j)] = ...
                lefatigue3(tau0,gamma0,X(i),dSL,Tf(j),tm,GIIth,GIIctm,G,gammaN,C1,mp,lambda(j),DeltaN,N,R,Lout,C(j),Ef,j,(Df/2),A(j));
        elseif i > nXcrit(j)
            leF(i,j,:)=Lout;
            lintact(i,j,:)=0;
            deltaA(i,j,:)=Lout;
            lcontrol(i,j,:)=Lout;
            iprop(i,j)=0;
            idebond(i,j)=0;
        end
    end
end
%=======Stress correspondent to the Gth value for each bundle level========

for i=1:n+1
    sigmath(i)=2*tau0/(lambda(i)*Tf(i))*sin(acos(1-gammath/gammaN));
end

for i=2:n+1
    leF(nXcrit(i)+1:end,i,1)=leF(nXcrit(i),i,1);
    
    lintact(nXcrit(i)+1:end,i,1)=lintact(nXcrit(i),i,1);
    deltaA(nXcrit(i)+1:end,i,1)=deltaA(nXcrit(i),i,1);
    lcontrol(nXcrit(i)+1:end,i,1)=lcontrol(nXcrit(i),i,1);
end
thindex=zeros(n+2,1);

for i=1:n+1
    for j=1:nX
        if leF(j,i,1)==leF(j,i,N)
            thindex(i)=j;
        else
            break
        end
    end
end

thindex(end)=thindex(end-1);

%=========================================================================

% % Bundle level 0 does not suffer from fatigue
% for i = 1:nX
%         leF(i,1,2:end)=leF(i,1,1);
%         lintact(i,1,2:end)=lintact(i,1,1);
%         deltaA(i,1,2:end)=deltaA(i,1,1);
%         lcontrol(i,1,2:end)=lcontrol(i,1,1);
% end



%=========================================================================

%=================Survival Probability (S.P) calculation===================
%Preallocating matrices
lnSur=zeros(nX,n+1,N);
lnSuc=zeros(nX,n+1,N);
lnSlr=zeros(nX,n+1,N);
lnSkr=zeros(nX,n+1,N);
lnSue=zeros(nX,n+1,N);
lnSue1=zeros(nX,n+1,N);
lnSke=zeros(nX,n+1,N);
lnSke1=zeros(nX,n+1,N);
lnSke2=zeros(nX,n+1,N);
% le=zeros(nX,n+1,N);
%-------------------------------------------------------------------------

% Infinite Fatigue Life definition;
% Runs the model for the specimen length, and uses the scaling law for the
% debonded phase
% leF(:,:,:)=Lout;
% lintact(:,:,:)=0;
% deltaA(:,:,:)=Lout;
% lcontrol(:,:,:)=Lout;

%============================S.P level-[0]=================================
for i=1:N
    if i==70
        1;
    end
    [lnSur(:,1,i),lnSue(:,1,i),lnSke1(:,1,i),lnSke2(:,1,i)]=...
        levelzerofinal(A(1),C(1),leF(:,1,i),lcontrol(:,1,i),deltaA(:,1,i),lintact(:,1,i));
end

%===========================S.P level-[1:n]================================

for i=1:n
    j=i+1;
    for h=1:N
        if j==2 && h==746
            1;
        end
        [lnSur(:,j,h),lnSuc(:,j,h)]=...
            levelifinalv2(lnSue(:,j-1,h),lnSke1(:,j-1,h),lnSke2(:,j-1,h),leF(:,j-1,h),lcontrol(:,j-1,h),deltaA(:,j-1,h),lintact(:,j-1,h));
    end
% Sur numerical control ===================================================
%Convert lnSur to a failure probability
    Suc(:,:)=exp(lnSuc(:,j,:));
    Suc(Suc>=1)=1-eps;

aux=Suc(:,N);
aux=1-aux;
erro1=1e-6;
if j==2
    for index=2:nX
        %     if abs(aux(index)-aux(index-1))<=erro
        %         indexLim=index;
        %     else
        %         break
        %     end
        if aux(index)<erro1
            indexLim=index;
        else
            break
        end
    end
end
        
%==========================================================================

    % Average of the caracteristic lengths for level j with the Failure probability associated to j, calculated previously
    leF1(:,:)=leF(:,j,:);
    lcontrol1(:,:)=lcontrol(:,j,:);
    lintact1(:,:)=lintact(:,j,:);
    deltaA1(:,:)=deltaA(:,j,:);
    
    %test the limits of this for loop, nXcrit or nX==========
    for stress=1:indexLim
        leF(stress,j,1:end)=leF1(stress,1);
        lcontrol(stress,j,1:end)=lcontrol1(stress,1);
        lintact(stress,j,1:end)=lintact1(stress,1);
        deltaA(stress,j,1:end)=deltaA1(stress,1);
    end
    
Suc(indexLim+1:end,:)=sort(Suc(indexLim+1:end,:),2,'descend');

erro2=30*eps;

    for stress=indexLim+1:nXcrit(j)
        aux1=[1, Suc(stress,1:end-1)];
        aux2=aux1-Suc(stress,:);
        aux2(aux2<erro2)=0;
%         aux2(1)=1-Suc(stress,1);
%         aux2=abs(aux2);
%         aux2=sort(aux2);
        aux3=conv(aux2,leF1(stress,:),'full');
        aux4=conv(aux2,lcontrol1(stress,:),'full');
        aux5=conv(aux2,lintact1(stress,:),'full');
        aux6=conv(aux2,deltaA1(stress,:),'full');
        aux31=aux3(1:N)./(1-Suc(stress,:));
        %
        aux31(aux31<=leF1(stress,1))=leF1(stress,1);
        leF(stress,j,:)=aux31;
        %
        aux41=aux4(1:N)./(1-Suc(stress,:));
        aux41(aux41<=lcontrol1(stress,1))=lcontrol1(stress,1);
        lcontrol(stress,j,:)=aux41;
        %
        aux51=aux5(1:N)./(1-Suc(stress,:));
        aux51(aux51<=lintact1(stress,1))=lintact1(stress,1);
        lintact(stress,j,:)=aux51;
        %
        aux61=aux6(1:N)./(1-Suc(stress,:));
        aux61(aux61<=deltaA1(stress,1))=deltaA1(stress,1);
        deltaA(stress,j,:)=aux61;
    end             

% test(:,:)=leF(:,j,:);
%==========================================================================
    for cycle=1:N
        lnSkrI1=zeros(nX,1);
        lnSkrI2=zeros(nX,1);
        lnSkeI1=zeros(nX,1);
        lnSkeI2=zeros(nX,1);
        lnSlrI=cumtrapz(lnSur(1:nX,j,cycle))*DX./X;
        lnSlrI(1)=0;
        if mod(k,1)==0
            %If k is integer, takes every other k-value:
            lnSkrI1(1:nK,1)=(k*lnSlrI(1:k:end)-lnSlrI(1:nK))./(k-1);
        else
            %If k is not integer (interpolation required):
            %Interpolates lnSlr for k*X (non-integer k):
            lnSlrKI(1:nK,1)=interp1(X,lnSlrI,XK,'linear');
            %uses lnSlr interpolated @ k*X (lnSlrK)
            lnSkrI1(1:nK,1)=(k*lnSlrKI-lnSlrI(1:nK))./(k-1);
        end
        
        lnSkrI1(nK+1:nX,1)=Ck*lnSur(nK+1:nX,j,cycle);
        lnSkrI1(1)=0;
        %----------------------------------------------------------------------
        %=========Survival probability calculated for the stress plateau=======
        lnSkrI2(1:nK)=lnSur(1:k:end,j,cycle);
        lnSkrI2(nK+1:end)=k^m*lnSur(nK+1:end,j,cycle);
        lnSkrI2(1)=0;
        %======================================================================
        lnSkrI=lnSkrI1+lnSkrI2;
        %%%For effective recovery length%%%%%%%%%%%%%%%%%%
        %Calculates survival functions
%         lnSke1(:,j,cycle)=lintact(:,j,cycle).*lnSkrI1;
%         lnSke2(:,j,cycle)=deltaA(:,j,cycle).*lnSkrI2;
        lnSke1(:,j,cycle)=leF(:,j,cycle).*lnSkrI1;
        lnSke2(:,j,cycle)=leF(:,j,cycle).*lnSkrI2;
        lnSue(:,j,cycle)=(leF(:,j,cycle)).*lnSur(:,j,cycle);
    end   
end

%========================Post-process results==============================
lnSuo=Lout/Lref*lnSur;
%Bundle strength distribution for output length--------------------------
Fuo=1-exp(lnSuo);
%-----------------Average bundle strength for output length----------------
intFuoDX=trapz(Fuo)*DX;
Xavgo=X(end)-intFuoDX;

%-----------------CoV of bundle strength for output length-----------------
intXFuoDX=zeros(1,n+1,N);
CoVo=sqrt(X(end).^2-Xavgo.^2-2.*intXFuoDX)./Xavgo;

%============Re-runs the model if free edges are to be considered==========
if and(T==41,freebounds==1)
    boundaries(n,A,C,Ccorner,Cedge);
end
end