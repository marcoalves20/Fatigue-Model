clc
clear 
n=18;
%d1=digits(50);
fatigueX(n,2)

global Fuo Xavgo X 
global N DeltaN
global sigmath nXcrit DX
global nX nK

%%
% figure()
% for i=1:n+1
%     plot(X,Fuo(:,i,1))
%     hold on
% %     legend('0','1', '2', '3', '4', '5', '6', '7')
% xlabel '\sigma^\infty (GPa)'
% ylabel 'F_U_r^[^i^] (%)'
% end
% 
% figure()
% plot([0:n],Xavgo(1,:,1))
% set(0,'defaultTextInterpreter','latex'); %trying to set the default
% xlabel 'level-[$i$]'
% ylabel '$X_m^{[i]}$'
% legend 'linear aprox.'
% aux1=[0:n];
% aux2=Xavgo(1,:,1);
% save('xdata.txt', 'aux1', '-ascii');
% save('ydata.txt', 'aux2', '-ascii');
%  
% % aux=Xavgo(1,:,1);
% % save('xmeziere.txt', 'aux', '-ascii');
% 
% figure()
% cycles=DeltaN*linspace(1,N,N);
% for i=1:N
% plot(X,Fuo(:,12,i))
% hold on
% end
% 
%  %% normal scale for different R values
% R1xdata=[1;27426.16034;33426.52321;69620.25316;1014767.932;1204641.35;1550632.911;1635021.097;508438.8186;259493.6709];
% R1ydata=[1; 0.75;0.65;0.600;0.497;0.495;0.49;0.487;0.510;0.543];
% % 

R1xdata=[1; 28605; 68595; 198505; 805485; 1005281];
R1ydata=[1; 0.7;0.62;0.573;0.5;0.463];

%Data to plot
staticexperimental=4916;
bundle=18;
cycles=zeros(length(DeltaN),1);
cycles(1)=1;
for i = 2:length(DeltaN)
    cycles(i)=cycles(i-1)+DeltaN(i);
end
Xaux(bundle,:)=Xavgo(1,bundle,1:end);
Strength(bundle,:)=Xaux(bundle,:);

% =====================plot template==========================
figure();
semilogx(R1xdata,staticexperimental*R1ydata,'kx')
hold on
Xaux(bundle,:)=Xavgo(1,bundle,1:end);
Strength(bundle,:)=Xaux(bundle,:);
% Strength(bundle,:)=Xaux(bundle,:);
% Strength(bundle,:)=Xaux(bundle,:);
cycles=zeros(length(DeltaN),1);
cycles(1)=DeltaN(1);

for i = 2:length(DeltaN)
    cycles(i)=cycles(i-1)+DeltaN(i);
end
% 
semilogx(cycles,Strength(bundle,:),'LineWidth',1.5)
hold on
% Set latex interpreter
set(0,'DefaultTextInterpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','fontsize',12)
legend ({'Experimental data', 'Model prediction'}, 'interpreter','latex')
legend boxoff
% ylim([0 1])
xlabel ('Number of cycles','fontsize',14)
ylabel ('$X_m^{[17]}$ (MPa)','fontsize',14)
% aux=Strength(bundle,:);
% aux2=cycles;
% save('12y.txt', 'aux', '-ascii');
% save('12x.txt', 'aux2', '-ascii');

%%
%Adaptive Step Plot

c1=[0,0,0]/255; % black
c2=[31,119,180]/255; %blueish
c3=[255,127,14]/255; % orange
c4=[44,160,44]/255; %green
figure()
set(0,'DefaultTextInterpreter', 'latex')
x0=20;
y0=15;
width=8;
height=7;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

Xaux(bundle,:)=Xavgo(1,bundle,1:end)/Xavgo(1,bundle,1);
Strength(bundle,:)=Xaux(bundle,:);
aux40=Strength(bundle,:);

save('aux40.txt', 'aux40', '-ascii');
save('cycles4.txt', 'cycles', '-ascii');

yyaxis left
semilogx(cycles,Strength(bundle,:),'Color',c2,'linewidth',1.5)
ylabel ('\% of Static Strength','fontsize',10)
set(gca,{'ycolor'},{'k'},'fontsize',11)
ylim([0.2 1])
xlim([1 100000000])
yyaxis right
semilogx(cycles,DeltaN,'--','Color',c2,'linewidth',2)
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
'Location','west');
legend boxoff
% 
% print('-painters','-depsc','adaptivestep')                  
%%
global leF deltaA
test(:,:)=deltaA(:,10,:);        
figure()
semilogx(cycles,test(35,:))
%% Save Files for crop
Xavgocrop(:,:)=Xavgo(1,:,:);
% save('cycles2.mat','cycles')
% save('Xavgo2.mat','Xavgocrop')
% save('leF2.mat','leF')
% save('deltaA2.mat','deltaA')