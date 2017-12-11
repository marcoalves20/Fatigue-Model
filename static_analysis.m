clc
clear all
nominalX(16,2)

global leF

%%
%effective recovery length comparison

global X
Ole=dlmread('myfile.txt');
lef=leF(:,:,1);
save('leF.txt', 'lef', '-ascii');
bundle=9;
Xindex=4921;

rgb_value1 = [14,34,228]/255;
rgb_value2 = [0,0,0]/255;

% Figure1=figure(1);clf;
% set(Figure1,'defaulttextinterpreter','latex');
figure()
plot(X(1:Xindex),Ole(1:Xindex,bundle),'Color',rgb_value1,'LineWidth',1)
hold on
plot(X(1:Xindex),lef(1:Xindex,bundle),'Color',rgb_value2,'LineWidth',1)
xlabel ('$\sigma^{\infty} \hspace{0.2mm} [MPa]$','fontsize',12)
ylabel ('$l_e/2 \hspace{0.2mm} [mm]$ ','fontsize',12)
legend('HFBM','Modified-HFBM','Location','northwest') 
legend boxoff  
title '$i=8$' 
%%
% avarage strength 

global Xavgo
Xavg(:)=Xavgo(:,:,1);
save('Xavg.txt', 'Xavg', '-ascii');
rgb_value1 = [14,34,228]/255;
rgb_value2 = [0,0,0]/255;

% Figure1=figure(1);clf;
% set(Figure1,'defaulttextinterpreter','latex');
figure()
plot([0:16],Xavg,'Color',rgb_value2,'LineWidth',1)
xlabel ('$level[i]$','fontsize',12)
ylabel ('$X_m^{[i]} (GPa)$ ','fontsize',12)




