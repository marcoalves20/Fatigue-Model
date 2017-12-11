% Life percentage in each phase
% clc 
% clear all
% n=17;
% fatigueX(n,2);
global iinit leF Lout
global Xavgo Xmax N
global DeltaN DX deltaA nX Lout
global sigmath
global lecrit
lecrit = 2*lecrit;

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
 l=nX;
 c=length(cyclesf);

maxstress=6800/DX;

%==========================================================================
%Determine the SN currve vector;
for i=1:n+1
    Xaux1(i,:)=Xavgof(i,1:end);
end

cycleRead=zeros(n+1,l);
for i = 1:n+1
    %for each stress value related to iinit
    for j=1:l
        auxvector=Xaux1(i,:);
        auxvector=auxvector(:);
        stress=j*DX;
        auxvector = auxvector - stress;
        [M,k] = min(abs(auxvector));
            if auxvector(k)==auxvector(length(cyclesf))
                k=length(cyclesf);
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
        if cycleRead(i,j)==0 | cycleRead(i,j) > 0.5e7
            cycleRead(i,j)=length(cyclesf);
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
         cycleRead(i,j)=cyclesf(aux);
        end
    end
end

% maxi=max(max(cycleRead));
for i =1:n+1
    for j =1:l
        if cycleRead(i,j)>5e6
            cycleRead(i,j)=Inf;
        end
    end
end
SN=cycleRead';
% SN(SN==max(max(SN)))=Inf
% SN(SN==1)=NaN;

%%
%==========================================================================
fenda(:,:)=leFf(:,2,:);
figure()
semilogx(cyclesf,fenda(18,:))

%%
%==========================================================================
% for i=1:nX
%     for j=1:n+1
%         for o=1:length(cyclesf)
%              if leFf(i,j,o)<=lecrit(j,1) %|| leFf(i,j,o)==Lout
%                 initiation(i,j)=cyclesf(o);
%              elseif leFf(i,j,o)>=Lout-1000*eps
%                  debonded(i,j)=cyclesf(o);
%                  continue
%             end
%         end
%     end
% end

for i=1:nX
    for j=1:n+1
        for o=1:length(cyclesf)
             if leFf(i,j,o)<=lecrit(i,j) %|| leFf(i,j,o)==Lout
                initiation(i,j)=cyclesf(o);
             elseif leFf(i,j,o)>=Lout-1000*eps
                 debonded(i,j)=cyclesf(o);
                 continue
            end
        end
    end
end


[t1,t2]=size(initiation);


for i=1:n+1   
    for j=1:t1
        initRatio(j,i)=initiation(j,i)/SN(j,18);
        debondRatio(j,i)=debonded(j,i)/SN(j,18);
    end
end

%Correction of the initiation Ratio

initRatio(initRatio>1)=1;
propRatio=1-initRatio;



        
% 
% [s1,s2]=size(SN);
% for j=s2
%     for i=1:s1
%         if initiation(i,j)>=SN(i,j);
%             initiation(i,j)=SN(i,j);
%         elseif initiation(i,j) < SN(i,j) 
%          propagation(i,j)=SN(i,j)-initiation(i,j);
%         end
%     end
% end  
% 
% initratio=initiation./SN;
% propratio=propagation./SN;


% debonded(debonded==Inf)=0;
%%
%==========================================================================
figure(1)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax1 = axes;
im1 = imagesc(log(initRatio(1:23,1:end)),'Parent',ax1);
set(ax1,'YDir','normal')
hcb=colorbar;
set(hcb,'YTick',[-7.36 -5.5225 -3.6850 -1.8475 -0.01],'YTickLabel',...
     {'0','25','50','75','100'})
ylabel(hcb,'\% of Fatigue Life','Interpreter', 'latex','fontsize',10)
hold on
aux1=17.5*ones(length([0:n+1]),1);
plot([0:n+1],aux1,'w--','linewidth',1.0)
text(17,20.5,'$\sigma_{th}$','Interpreter', 'latex','fontsize',11,'Color','w')
% 'YTick',[1,3.98,6.99,8.75]
colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16 19],'Xticklabels',{'0', '3' '6', '9', '12', '15','18'},'YTick',[0 10 20 30 40 50],'Yticklabels',{'0', '1' '2', '3', '4', '5'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
% hold on
% plot([1:n+1],ydata(2:n+2)/DX,'Color',Color1,'linewidth',1.2)


%% 
%==========================================================================
figure(2)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax2 = axes;
im2 = imagesc(log(propRatio(1:23,1:end)),'Parent',ax2);
set(ax2,'YDir','normal')
% hcb=colorbar('YTickLabel',...
%     {'0','0.2','0.4','0.6','0.8','1'});
hcb=colorbar();
% set(hcb,'TickLabelInterpreter','latex','YTickLabel',...
%      {'0','0.2','0.4','0.6','0.8','1'})
set(hcb,'YTick',[-3 -2.25 -1.5 -0.75 -0.01],'YTickLabel',...
     {'0','25','50','75','100'})
ylabel(hcb,'\% of Fatigue Life','Interpreter', 'latex','fontsize',10)
hold on
% plot([0:n+1],aux1,'w--','linewidth',1.0)
aux1=17.5*ones(length([0:n+1]),1);
plot([0:n+1],aux1,'w--','linewidth',1.0)
text(17,20.5,'$\sigma_{th}$','Interpreter', 'latex','fontsize',11,'Color','w')
colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16 19],'Xticklabels',{'0', '3' '6', '9', '12', '15','18'},'YTick',[0 10 20 30 40 50],'Yticklabels',{'0', '1' '2', '3', '4', '5'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
% hold on
% plot([1:n+1],ydata(2:n+2)/DX,'Color',Color1,'linewidth',1.2)

%%
figure(3)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax3 = axes;
im3 = imagesc(log(debondRatio(1:49,2:end)),'Parent',ax3);
set(ax3,'YDir','normal')
% hcb=colorbar('YTickLabel',...
%     {'0','0.2','0.4','0.6','0.8','1'});
hcb=colorbar();
% set(hcb,'TickLabelInterpreter','latex','YTickLabel',...
%      {'0','0.2','0.4','0.6','0.8','1'})
set(hcb,'YTick',[-3 -2.25 -1.5 -0.75 -0.01],'YTickLabel',...
     {'0','0.25','0.5','0.75','1'})
ylabel(hcb,'\% of Fatigue Life','Interpreter', 'latex','fontsize',10)

colormap(cmap);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[0 3 6 9 12 15 18],'YTick',[0 10 20 30 40 50],'Yticklabels',{'0', '1' '2', '3', '4', '5'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
%==========================================================================
%%
figure(7);
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax7 = axes;
im7 = imagesc(log(SN(1:35,2:end)),'Parent',ax7);
set(ax7,'YDir','normal')
hcb=colorbar();
set(hcb,'YTick',[0.75 4.375 8 11.625 15.25],'YTickLabel',...
     {'0','4','8','12','$\infty$'})
%'YTickLabel',{'0','6','8','10','12','14','$\infty$'}
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'$\log$(N. cycles)','Interpreter', 'latex','fontsize',10)
colormap(cmap);
colormap(flipud(colormap))
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[1 4 7 10 13 16],'XTicklabels',[1 4 7 10 13 16],'YTick',[0 10 20 30 40 50 60 70],'Yticklabels',{'0', '1' '2', '3', '4', '5','6', '7'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'
hold on 
plot([0:n],Xavgo(1,:,1)/DX,'k','linewidth',1.2)
% hold on
% plot([0:n],sigmath,'Color',Color1,'linewidth',1.2)
% hold on
plot([0:n],ydata(1:n+1)/DX,'Color',Color1,'linewidth',1.2)
hold on
aux1=17.5*ones(length([0:n+1]),1);
plot([0:n+1],aux1,'w--','linewidth',1.0)
legend({'$X_m^{[i]}$','$X_{fd}^{[i]}$','$\sigma_{th}$'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'Color',[1 1 1],...
'Location','east');
% legend boxoff ('$\sigma^{\infty}_{crt}$,'interpreter','latex')
legend boxoff
% 
%             
