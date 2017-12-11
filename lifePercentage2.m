% Life percentage in each phase
clc 
clear all
n=17;
fatigueX(n,2);

global iinit idebond iprop
global Xavgo
global DeltaN DX
global sigmath
sigmath=sigmath/10;
%%
global nXcrit
criticalstress=nXcrit;
Color1=[255 127 14]/255;
ydata=dlmread('ydata.txt');
%==========================================================================
%My custom colormap
c1=[48 58 147]/255; %Blue
c2=[0 175 235]/255; %Light blue
c3=[38 176 76]/255; %Green
c4=[250 239 20]/255; % Yellow
c5=[240 79 39]/255; %Red
n1=20;

cmap=[linspace(c1(1),c2(1),n1);linspace(c1(2),c2(2),n1);linspace(c1(3),c2(3),n1)];
cmap(:,n1+1:2*n1)=[linspace(c2(1),c3(1),n1);linspace(c2(2),c3(2),n1);linspace(c2(3),c3(3),n1)];
cmap(:,2*n1+1:3*n1)=[linspace(c3(1),c4(1),n1);linspace(c3(2),c4(2),n1);linspace(c3(3),c4(3),n1)];
cmap(:,3*n1+1:4*n1)=[linspace(c4(1),c5(1),n1);linspace(c4(2),c5(2),n1);linspace(c4(3),c5(3),n1)];
cmap=cmap';
%==========================================================================
[l,c]=size(iinit);

maxstress=6800/DX;
cycles=zeros(length(DeltaN),1);
cycles(1)=1;
for i = 2:length(DeltaN)
    cycles(i)=cycles(i-1)+DeltaN(i);
end

%==========================================================================
%==========================================================================

%==========================================================================
%Determine the SN currve vector;
for i=1:n+1
    Xaux(i,:)=Xavgo(1,i,1:end);
end

cycleRead=zeros(n+1,l);
for i = 1:n+1
    %for each stress value related to iinit
    for j=1:l
        auxvector=Xaux(i,:);
        auxvector=auxvector(:);
        stress=j*DX;
        auxvector = auxvector - stress;
        [M,k] = min(abs(auxvector));
            if auxvector(k)==auxvector(length(DeltaN))
                k=length(DeltaN);
%             elseif i_row==0
%                     cycleRead(i,j)=1;
            else
                [i_row, i_col] = ind2sub(size(auxvector),k);
                cycleRead(i,j)=i_row;
           end
    end
end

for i=1:n+1
    for j=1:l
        if cycleRead(i,j)==0 | cycleRead(i,j) > 1e7
            cycleRead(i,j)=length(DeltaN);
        end
    end
end

%Convert index of cycleRead to actual number of cycles
for i=1:n+1
    for j=1:l
        if cycleRead(i,j)==0;
            1;
        else
            aux=cycleRead(i,j);
         cycleRead(i,j)=cycles(aux);
        end
    end
end

maxi=max(max(cycleRead));
for i =1:n+1
    for j =1:l
        if cycleRead(i,j)>1e7
            cycleRead(i,j)=Inf;
        end
    end
end
SN=cycleRead';
%==========================================================================
debonded=zeros(l,n+1);
propagation=zeros(l,n+1);


for i =1:l
    for j=1:n+1
        if iinit(i,j) >= SN(i,j)
            iinit(i,j)=SN(i,j);
            iprop(i,j)=0;
            idebond(i,j)=0;
        else
            if iinit(i,j)+iprop(i,j)>=SN(i,j)
                iprop(i,j)=SN(i,j)-iinit(i,j);
                idebond(i,j)=0;
            end
            if iinit(i,j)+iprop(i,j) <= SN(i,j) 
                idebond(i,j)=SN(i,j)-iinit(i,j)-iprop(i,j);
            end
        end            
    end
end


% debonded(debonded==Inf)=0;

%==========================================================================
figure(1)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax1 = axes;
im1 = imagesc(log(iinit(1:maxstress,2:end)),'Parent',ax1);
set(ax1,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex','YTick',[10 12 14 16])
ylabel(hcb,'log of N. Cycles','Interpreter', 'latex','fontsize',10)
% 'YTick',[1,3.98,6.99,8.75]
colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500 600],'Yticklabels',{'0', '1' '2', '3', '4', '5', '6'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
hold on
plot([1:n+1],ydata(2:n+2)/DX,'Color',Color1,'linewidth',1.2)


%==========================================================================
figure(2)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax2 = axes;
im2 = imagesc(log(iprop(1:maxstress,2:end)),'Parent',ax2);
set(ax2,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex','YTick',[10 12 14 16])
ylabel(hcb,'log of N. Cycles','Interpreter', 'latex','fontsize',10)
% 'YTick',[1,3.98,6.99,8.75]
colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16 19],'XTicklabels',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500 600],'Yticklabels',{'0', '1' '2', '3', '4', '5','6'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
hold on
plot([1:n+1],ydata(2:n+2)/DX,'Color',Color1,'linewidth',1.2)

%==========================================================================
figure(3)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])
% 
% ratioplot=transpose(lifeRatio(1:end,1:maxstress));
ax3 = axes;
im3 = imagesc(log(idebond(1:maxstress,:)),'Parent',ax3);
set(ax3,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'$\%$ of Fatigue Life','Interpreter', 'latex','fontsize',10)
colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16 19],'XTicklabels',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500 600],'Yticklabels',{'0', '1' '2', '3', '4', '5','6'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
hold on 
plot([1:n+1],Xavgo(1,:,1)/DX)
hold on
plot([1:n+1],ydata(2:n+2)/DX,'Color',Color1,'linewidth',1.2)
%%
%==========================================================================

initiationRatio=iinit./SN;
figure(4)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])
% 
% ratioplot=transpose(lifeRatio(1:end,1:maxstress));
ax4 = axes;
im4 = imagesc(initiationRatio(1:maxstress,:),'Parent',ax4);
set(ax4,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'$\%$ of Fatigue Life','Interpreter', 'latex','fontsize',10)
colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16 19],'XTicklabels',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500 600],'Yticklabels',{'0', '1' '2', '3', '4', '5','6'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
hold on
plot([0:n],Xavgo(1,:,1)/DX,'k','linewidth',1.2)
hold on
plot([1:n+1],criticalstress,'w','linewidth',1.2)
hold on
plot([1:n+1],ydata(2:n+2)/DX,'Color',Color1,'linewidth',1.2)
legend({'$X_m^{[i]}$','$\sigma^{\infty}_{crt}$','$\sigma^{\infty}_{th}$'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'TextColor', [1 1 1],...
'Location','east');
% legend boxoff ('$\sigma^{\infty}_{crt}$,'interpreter','latex')
legend boxoff
%==========================================================================


propRatio=iprop./SN;
figure(5)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])
% 
% ratioplot=transpose(lifeRatio(1:end,1:maxstress));
ax5 = axes;
im5 = imagesc(propRatio(1:maxstress,:),'Parent',ax5);
set(ax5,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'$\%$ of Fatigue Life','Interpreter', 'latex','fontsize',10)
colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16 19],'XTicklabels',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500 600],'Yticklabels',{'0', '1' '2', '3', '4', '5','6'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
hold on 
plot([0:n],Xavgo(1,:,1)/DX,'k','linewidth',1.2)
hold on
plot([1:n+1],criticalstress,'w','linewidth',1.2)
hold on
plot([0:n],ydata(1:n+1)/DX,'Color',Color1,'linewidth',1.2)
legend({'$X_m^{[i]}$','$\sigma^{\infty}_{crt}$','$\sigma^{\infty}_{th}$'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'TextColor', [1 1 1],...
'Location','east');
% legend boxoff ('$\sigma^{\infty}_{crt}$,'interpreter','latex')
legend boxoff

%==========================================================================

debondRatio=idebond./SN;
figure(6)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])
% 
% ratioplot=transpose(lifeRatio(1:end,1:maxstress));
ax6 = axes;
im6 = imagesc(debondRatio(1:maxstress,:),'Parent',ax6);
set(ax6,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'$\%$ of Fatigue Life','Interpreter', 'latex','fontsize',10)
colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16 19],'XTicklabels',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500 600],'Yticklabels',{'0', '1' '2', '3', '4', '5','6'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
hold on
plot([0:n],Xavgo(1,:,1)/DX,'k','linewidth',1.2)
hold on 
plot([1:n+1],criticalstress,'w','linewidth',1.2)
hold on
plot([0:n],ydata(1:n+1)/DX,'Color',Color1,'linewidth',1.2)
legend({'$X_m^{[i]}$','$\sigma^{\infty}_{crt}$','$\sigma^{\infty}_{th}$'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'Color',[1 1 1],...
'Location','east');
% legend boxoff ('$\sigma^{\infty}_{crt}$,'interpreter','latex')
legend boxoff

% hold on
% plot([1:n+1],sigmath,'','linewidth',1.2)

%==========================================================================

figure(7);
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax7 = axes;
im7 = imagesc(log(SN(1:maxstress,2:end)),'Parent',ax7);
set(ax7,'YDir','normal')
hcb=colorbar('YTickLabel',...
    {'0','10','12','$\infty$',...
     '$\infty$'});
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'log(N. of cycles)','Interpreter', 'latex','fontsize',10)
colormap(cmap);
colormap(flipud(colormap))
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16 19],'XTicklabels',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500 600],'Yticklabels',{'0', '1' '2', '3', '4', '5','6'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
hold on 
plot([0:n],Xavgo(1,:,1)/DX,'k','linewidth',1.2)
hold on
% plot([0:n],sigmath,'Color',Color1,'linewidth',1.2)
% hold on
plot([0:n],ydata(1:n+1)/DX,'Color',Color1,'linewidth',1.2)
% hold on
% plot([1:n+1],criticalstress,'w','linewidth',1.2)
legend({'$X_m^{[i]}$','$\sigma^{\infty}_{th}$'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'Color',[1 1 1],...
'Location','east');
% legend boxoff ('$\sigma^{\infty}_{crt}$,'interpreter','latex')
legend boxoff

            
