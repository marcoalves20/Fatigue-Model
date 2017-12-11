% Attempt to plot a fatigue life diagram, by splitting the run of the model
% in 2. Obtain 2 different matrices, and combine them togeather. Fuo1 and
% Fuo2

Fuo1=dlmread('Fuo1.txt');
Fuo2=dlmread('Fuo2.txt');
Fuo3=dlmread('Fuo3.txt');
%%
global DX
Fuofinal=zeros(4001,10000+2323+20-23);
Fuofinal(:,1:20)=Fuo1(:,:);
Fuofinal(:,21:2343)=Fuo2(:,2:end);
Fuofinal(:,2344:end)=Fuo3(:,24:end);
intFuoDX=trapz(Fuofinal)*DX;
Xavgo=X(end)-intFuoDX;
% 
cycles=zeros(length(Xavgo),1);
cycles(1:20)=linspace(1,20,20);

for i = 21:2343
    cycles(i)=cycles(i-1)+10;
end
for i = 2344:length(Xavgo)
    cycles(i)=cycles(i-1)+1000;
end
% 
% % 
figure()
semilogx(cycles, Xavgo)
%%
probMatrix(:,:)= Fuofinal(1:381,1:end);
x=linspace(1,10e7,12320);
y=linspace(0.2,1,381);
[cs,hc]=contourf(x,y,probMatrix);

set(hc,...
    'EdgeColor','none')
set(gca,'xscale','log')