% 1st part of the SN curve 
clc

load('cycles1.mat')
cycles1=cycles;
load('Xavgo1.mat')
Xavgo1=Xavgocrop;
load('leF1.mat')
leF1=leF;
load('deltaA1.mat')
deltaA1=deltaA;

%%
%2nd part of the SN curve
load('cycles2.mat')
cycles2=cycles;
load('Xavgo2.mat')
Xavgo2=Xavgocrop;
load('leF2.mat')
leF2=leF;
load('deltaA2.mat')
deltaA2=deltaA;
%%
x2=301;
aux=cycles2(x2:end);

% Combine SN curve
cyclesf=zeros(length(cycles1)+length(cycles2)-x2+1,1);
cyclesf(1:length(cycles1))=cycles1(:,1);
cyclesf(length(cycles1)+1:end)=cycles2(x2:end,1);

Xavgof=zeros(19,length(cycles1)+length(cycles2)-x2+1);
Xavgof(:,1:length(cycles1))=Xavgo1(:,:);
Xavgof(:,length(cycles1)+1:end)=Xavgo2(:,x2:end);

leFf=zeros(nX,19,length(cycles1)+length(cycles2)-x2+1);
leFf(:,:,1:length(cycles1))=leF1(:,:,:);
leFf(:,:,length(cycles1)+1:end)=leF2(:,:,x2:end);

deltaAf=zeros(nX,19,length(cycles1)+length(cycles2)-x2+1);
deltaAf(:,:,1:length(cycles1))=deltaA1(:,:,:);
deltaAf(:,:,length(cycles1)+1:end)=deltaA2(:,:,x2:end);

%% 
%Full SN curve
figure()
semilogx(cyclesf,Xavgof(18,:))
auxx=Xavgof(18,:);
save('valx.txt', 'cyclesf', '-ascii');
save('valy.txt', 'auxx', '-ascii');
%%
%Adaptive step
figure()
set(0,'DefaultTextInterpreter', 'latex')
x0=20;
y0=15;
width=8;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
DeltaN=zeros(length(cyclesf),1);
DeltaN(1:1000)=1;
DeltaN(1001:32425)=100;
DeltaN(32426:end)=1000;
yyaxis left
semilogx(cyclesf,Xavgof(18,:)/Xavgof(18,1)*100,'Color',c2,'linewidth',1.5)
ylabel ('\% of Static Strength','fontsize',10)
set(gca,{'ycolor'},{'k'},'fontsize',11)
ylim([20 100])
xlim([1 100000000])
yyaxis right
semilogx(cyclesf,DeltaN,'--','Color',c3,'linewidth',2)
ylabel ('$\Delta N$','fontsize',10)
rax = findall(0, 'YAxisLocation','right');
lax = findall(0, 'YAxisLocation','left');
set(lax,'Yscale','linear','box','off');
set(rax,'Yscale','log','box','off');
xlabel ('Number of cycles','fontsize',11)

set(gca,'TickLabelInterpreter','latex','FontSize',10,'Xtick',[1 100 10000 1000000 100000000],'Ytick',[1  10 100 1000 10000])
set(gca,{'ycolor'},{'k'},'fontsize',11)  % Left color red, right color blue...
set(gca, 'YLim', [1 10000])
h =legend({'Model','$\Delta N$ '},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',10,...
'Position',[0.26 0.6 0.1 0.2]);
legend boxoff

%%
test2(:,:)=leFf(:,2,:);        
figure()
semilogx(cyclesf,test2(20,:))




