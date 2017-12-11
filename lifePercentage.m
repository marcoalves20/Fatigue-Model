% Life percentage in each phase
clc 
clear all
n=16;
fatigueX(16,2);
%%
global iprop idebond
global Xavgo
global DeltaN DX
diff=idebond-iprop;
%%
% -------------------------------------------------------------------------------------------------------------------------------------------------
maxstress=500;
cycles=zeros(length(DeltaN),1);
cycles(1)=1;
for i = 2:length(DeltaN)
    cycles(i)=cycles(i-1)+DeltaN(i);
end

[l,c]=size(iprop);
% for i=1:l
%     for j=1:c
%         if iprop(i,j)~=0
%             iprop(i,j)=cycles(iprop(i,j),1);
%         end
%     end
% end

%------------------------
%Figure 1 - duration of the crack initiation phase, no of cycles
figure(1)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax1 = axes;
im1 = imagesc(log(iprop(1:maxstress,2:end)),'Parent',ax1);
set(ax1,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex','YTick',[10 12 14 16])
ylabel(hcb,'log of N. Cycles','Interpreter', 'latex','fontsize',10)
% 'YTick',[1,3.98,6.99,8.75]
colormap(jet);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500],'Yticklabels',{'0', '1' '2', '3', '4', '5'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'


%==========================================================================

% plot the percentage of the fatigue life in each phase;
%1st have to get the SN curve, and for each stress value, check the number
%of cycles to failure and compare to the correspondent value in iprop, and
%idebond.
for i=1:n+1
    Xaux(i,:)=Xavgo(1,i,1:end);
end

cycleRead=zeros(n+1,l);
for i = 1:n+1
    %for each stress value related to iprop
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
            cycleRead(i,j)=maxi;
        end
    end
end
SN=cycleRead';

%Get lifeRatio between iprop and total fatigue life (iprop/cycleRead)
lifeRatio=zeros(n+1,l);
for i =1:n+1
    for j=1:l 
        lifeRatio(i,j)=iprop(j,i)./cycleRead(i,j);
        if lifeRatio(i,j) > 1
            lifeRatio(i,j) = 1;
        end
    end
end

% Fig 3 - % of initiation phase in total fatigue life =====================
% figure(3)
% x0=20;
% y0=15;
% width=9.5;
% height=7;
% set(0,'DefaultTextInterpreter', 'latex')
% set(gcf,'units','centimeters','position',[x0,y0,width,height])
% 
% ratioplot=transpose(lifeRatio(1:end,1:maxstress));
% ax3 = axes;
% im3 = imagesc(ratioplot,'Parent',ax3);
% set(ax3,'YDir','normal')
% hcb=colorbar;
% set(hcb,'TickLabelInterpreter','latex')
% ylabel(hcb,'$\%$ of Fatigue Life','Interpreter', 'latex','fontsize',10)
% colormap(jet);
% set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500],'Yticklabels',{'0', '1' '2', '3', '4', '5'});
% ylabel '$\sigma^{\infty}$ (GPa)'
% xlabel '$i$'

%2 version
initiationRatio=iprop./SN;
for i = 1:l
    for j=1:n+1
        if initiationRatio(i,j)>1
            initiationRatio(i,j)=1;
        end
    end
end

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
im3 = imagesc(initiationRatio(1:maxstress,:),'Parent',ax3);
set(ax3,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'$\%$ of Fatigue Life','Interpreter', 'latex','fontsize',10)
colormap(jet);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500],'Yticklabels',{'0', '1' '2', '3', '4', '5'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'


%==========================================================================
% Plot percentage of the propagation phase
% for i=1:l
%     for j=1:c
%         if idebond(i,j)~=0
%             idebond(i,j)=cycles(idebond(i,j),1);
%         end
%     end
% end
% 
% for i =1:1:l
%     for j =n+1
%         idebond(i,j)=idebond(i,j)-iprop(i,j);
% %         if idebond(i,j) < 0
% %             idebond(i,j)=0;
% %         end
%     end
% end


for i=1:l
    for j=1:c
         if diff(i,j) < 0
           diff(i,j)=0;
         elseif diff(i,j)~=0
            diff(i,j)=cycles(diff(i,j),1);
        end
    end
end

debondRatio=zeros(n+1,l);
for i =1:n+1
    for j=1:l 
        debondRatio(i,j)=diff(j,i)./cycleRead(i,j);
        if debondRatio(i,j) > 1
            debondRatio(i,j) = 1;
        end
    end
end
%%
% figure(4)
% x0=20;
% y0=15;
% width=9.5;
% height=7;
% set(0,'DefaultTextInterpreter', 'latex')
% set(gcf,'units','centimeters','position',[x0,y0,width,height])
% 
% ratioplot=transpose(debondRatio(2:end,1:maxstress));
% ax3 = axes;
% im3 = imagesc(ratioplot,'Parent',ax3);
% set(ax3,'YDir','normal')
% hcb=colorbar;
% set(hcb,'TickLabelInterpreter','latex')
% ylabel(hcb,'$\%$ of Fatigue Life','Interpreter', 'latex','fontsize',10)
% colormap(jet);
% set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500],'Yticklabels',{'0', '1' '2', '3', '4', '5'});
% ylabel '$\sigma^{\infty}$ (GPa)'
% xlabel '$i$'

figure(4)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ratioplot=transpose(debondRatio(2:end,1:maxstress));
ax3 = axes;
im3 = imagesc(ratioplot,'Parent',ax3);
set(ax3,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'$\%$ of Fatigue Life','Interpreter', 'latex','fontsize',10)
colormap(jet);
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[0 3 6 9 12 15 18],'YTick',[0 100 200 300 400 500],'Yticklabels',{'0', '1' '2', '3', '4', '5'});
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'

%Absolute values=========================================================
figure(5)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

transdiff=(idebond(1:maxstress,:));
ax5 = axes;
im5 = imagesc(log(idebond(1:maxstress,:)),'Parent',ax5);
set(ax5,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'No of cycles','Interpreter', 'latex','fontsize',10)
colormap(jet);
set(gca,'TickLabelInterpreter','latex','fontsize',10);
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'

%==========================================================================
%SN curve in a colormap plot


test=cycleRead';
figure(6)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])

ax6 = axes;
im6 = imagesc(log(test(1:maxstress,:)),'Parent',ax6);
set(ax6,'YDir','normal')
hcb=colorbar;
set(hcb,'TickLabelInterpreter','latex')
ylabel(hcb,'No of cycles','Interpreter', 'latex','fontsize',10)
colormap(jet);
set(gca,'TickLabelInterpreter','latex','fontsize',10);
ylabel '$\sigma^{\infty}$ (GPa)'
xlabel '$i$'


%==========================================================================

%==========================================================================
%3D plot of the S-N curve

figure(7)
x0=20;
y0=15;
width=9.5;
height=7;
set(0,'DefaultTextInterpreter', 'latex')
set(gcf,'units','centimeters','position',[x0,y0,width,height])
ax7 = axes;
surf(cycles,[1:16],Xaux(2:end,:),'EdgeColor','none')
set(ax7,'XScale','log')
xlabel ('N. of cycles','fontsize',11)
ylabel ('$i$','fontsize',11)
zlabel('$X_m^{[i]}$ (MPa)')%==========================================================================
set(gca,'TickLabelInterpreter','latex','fontsize',10,'XTick',[10 1000 100000 10000000],'YTick',[0 3 6 9 12 15 18]);
zlim([1500 7000])
ylim([1 18])


        
