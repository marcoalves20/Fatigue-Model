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
%3rd part of the SN curve
load('cycles3.mat')
cycles3=cycles;
load('Xavgo3.mat')
Xavgo3=Xavgocrop;
load('leF2.mat')
leF3=leF;
load('deltaA3.mat')
deltaA3=deltaA;

%% Crop paramters
x2=101;
x3=201
aux=cycles2(x2:end);

%%
% Combine SN curve
cyclesf=zeros(length(cycles1)+length(cycles2)-x2+1+length(cycles3)-x3+1,1);
cyclesf(1:length(cycles1))=cycles1(:,1);
cyclesf(length(cycles1)+1:length(cycles1)+length(cycles2(x2:end,1)))=cycles2(x2:end,1);
cyclesf(length(cycles1)+length(cycles2(x2:end,1))+1:end)=cycles3(x3:end,1);


Xavgof=zeros(n+1,length(cycles1)+length(cycles2)-x2+1+length(cycles3)-x3+1);
Xavgof(:,1:length(cycles1))=Xavgo1(:,:);
Xavgof(:,length(cycles1)+1:length(cycles1)+length(cycles2(x2:end,1)))=Xavgo2(:,x2:end);
Xavgof(:,length(cycles1)+length(cycles2(x2:end,1))+1:end)=Xavgo3(:,x3:end);
bundle=18;
figure()
semilogx(cyclesf,Xavgof(bundle,:))

aux=cyclesf;
aux2=Xavgof(bundle,:);
%%
save('comb4x.txt', 'aux', '-ascii');
save('comb4y.txt', 'aux2', '-ascii');

