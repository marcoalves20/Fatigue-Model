clc
clear all
n=18;
fatigueX(n,2);
global Fuo N DeltaN nX X Xavgo DX
%Choose the bundle to calculate the SN curve
bundle=17;
%Probability associated to the SN curve
probability=0.5;
%Static strength for the chosen bundle
Xstatic=Xavgo(1,bundle,1);
%%
% % for i=1:18
% bundle=10
% Figure1=figure(1);clf;
% set(Figure1,'defaulttextinterpreter','latex');
% pr=[0;0.1;0.2;0.3:0.4;0.5;0.6;0.7;0.8; 0.9;0.99];
% for k=1:length(pr)
% for i=1:nX
%         x=Fuo(i,bundle,:)-pr(k);
%         y=[1:N];
%         point=0;
%         [x, index] = unique(x); 
%         if size(x)==[1 1]
%             x(2)=x(1);
%             index(2)=index(1);
%         else
%         Nbreak(i)=(interp1(x,y(index),point))*DeltaN;
%         end
% end
% dim=size(Nbreak);
% Xaux=transpose(X(1:dim(2)));
% plot(Nbreak,Xaux,'LineWidth',1);
% hold on
% % Matrix(k,:,:)
% aux1(:,k)=Nbreak;
% aux2(:,k)=Xaux;
% end
% 
% R1xdata=[10799.225; 12657.14; 15045.04; 21974.99; 33538.56; 51537.25; 35111.92; 50352.4;73454.85];
% R1ydata=[0.801;0.799;0.799;0.7;0.7;0.7;0.67;0.669;0.67];
% staticexperimental=3500;
% expData=dlmread('expData.txt');
% dataX=expData(:,1);
% dataY=staticexperimental*expData(:,2);
% plot(R1xdata,staticexperimental*R1ydata,'kx')
% xlabel 'Number of cycles'
% ylabel '$X_m$'
% xlim([0 80000])
% ylim([2100 3500])
% legend('Failure Probability = 10%','Failure Probability = 50%','Failure Probability = 90%')
% legend boxoff
%% Probability color map
% bundle = 18;
% DeltaN;
% maxstress=410;
% minstress=1;
% probabilityMatrix(:,:) = Fuo(minstress:maxstress,bundle,:);
% aux=size(probabilityMatrix);
% x = [1 10^5];
% y = [0 1];
% 
% % plotting
% % Figure1=figure(1);clf;
% % set(Figure1,'defaulttextinterpreter','latex');
% R1Xdata=[27426.16034;29426.52321;69620.25316;1014767.932;1204641.35;1550632.911;1635021.097;508438.8186;259493.6709];
% R1ydata=[0.75;0.65;0.600;0.497;0.494;0.48;0.478;0.516;0.543];
% % hold on
% image(probabilityMatrix)
% image(x,y,probabilityMatrix,'CDataMapping','scaled')
% set(gca,'YDir','normal')
% colorbar
% colormap(parula)
% hold on
% plot(R1xdata,R1ydata,'mo',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','k',...
%     'MarkerSize',2)
% hold on
% Xaux(bundle,:)=Xavgo(1,bundle,1:end);
% Strength(bundle,:)=Xaux(bundle,:)/(Xavgo(1,bundle,1));
% cycles=linspace(1,aux(2),aux(2));
% plot(DeltaN*cycles(1,2:end),Strength(bundle,2:end),'k','LineWidth',1.5)
% 
% h = legend('Experimental Data','Model prediction');
% set(h,'color','none')
% legend boxoff
% xlabel('Number of Cycles','Interpreter','latex');
% ylabel('Strength (MPa)','Interpreter','latex');
%%
%Contour plot
maxstress=55;
minstress=25;
probabilityMatrix(:,:) = Fuo(minstress:maxstress,bundle,:);
staticexperimental=4916;

[l,c] = size(probabilityMatrix);
test=zeros(l,c);
for i = 1:l
    for j = 1:c
        if probabilityMatrix(i,j) < 0.001
            probabilityMatrix(i,j)=0;
        end
    end
end

cycles=zeros(length(DeltaN),1);
cycles(1)=DeltaN(1);
for i = 2:length(DeltaN)
    cycles(i)=cycles(i-1)+DeltaN(i);
end
Xaux(bundle,:)=Xavgo(1,bundle,1:end);
Strength(bundle,:)=Xaux(bundle,:);
Yvalue=linspace(2350,5500,l);

[hC hC] = contourf(cycles,Yvalue,probabilityMatrix);
set(hC,'LineStyle','none');
set(gca,'xscale','log')
colormap();
colorbar;
hold on
% R1Xdata=[27426.16034;29426.52321;69620.25316;1014767.932;1204641.35;1550632.911;1635021.097;508438.8186;259493.6709];
% R1ydata=[0.75;0.65;0.600;0.497;0.494;0.48;0.478;0.516;0.543];

R1Xdata=[ 26605; 65595; 198505; 805485; 1005281];
R1ydata=[0.7;0.62;0.573;0.5;0.463];

semilogx(R1Xdata,staticexperimental*R1ydata,'kx')
hold on
aux(bundle,:)=Xavgo(1,bundle,1:end);
Strength(bundle,:)=Xaux(bundle,:);
semilogx(cycles,Strength(bundle,:),'k','LineWidth',1.5)


