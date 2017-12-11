clc
clear all
global lnSuo Lout Lref lnSur N DeltaN DX X ER dAdN nX GIIctm GIIth

n=15;
fatigueX(n,2);

bundle=1;
nx=25;
Survival=zeros(N,1);

Survival(:)=(lnSuo(nx,bundle,:));

% plot([1:N],Survival)

%%
% n=7;
% fatigueX(n,2);
% 
% lnSuo=Lout/Lref*lnSur;
% 
% sigmamax=2000;
% nXc=sigmamax/DX+1;
% 
% FP=1-exp(lnSuo);
% Failure=FP(1:nXc,:,:);
% 
% intFuoDX=trapz(Failure)*DX;
% Xavgo=X(nXc)-intFuoDX;
% 
% bundle=5;
% Xbav(1:N)=Xavgo(1,bundle,:);
% 
% figure()
% plot(Xbav/Xbav(1))

%%
%Plot energy release rate vs crack growth rate

% dadn=zeros(nX,1);
% G=zeros(nX,1);

dadn(:)=dAdN(8,15,:);
G(:)=ER(8,15,:)/GIIctm;

plot(G,dadn)
xlabel 'G/G_c'
ylabel 'da/dN [mm/cycle]' 

%%
% WLT
n=15;
WLT(n);
global FuoWLT XavgoWLT
global Xavgo
figure()
plot([0:n],XavgoWLT(1,:))




