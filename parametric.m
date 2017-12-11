function[] = parametric(n)

global variAll variFew
global statsAll statsFew XFew
global lnSuoFew0 lnSuoFew1 lnSuoFew2 lnSuoFew5 lnSuoFew10 lnSuoFew15 lnSuoFew20  
global Xavgo CoVo lnSuo X nX
global varinow ivari
global nXFew iFew Xjumpfac Xmaxfac


%variAllTsl=1:1:150; 
%variFewTsl=[1,5,10,20,30,50,70,100,150]; 
%variAllTslAsym=[(0.0001:0.0001:0.001),(0.002:0.001:0.01),(0.02:0.01:0.1),(0.2:0.1:1),(2:1:10),(20:10:100),(200:100:1000),(2000:1000:10000),(20000:10000:100000),(200000:100000:1000000),(2000000:1000000:10000000),(20000000:10000000:100000000),(200000000:100000000:1000000000),(2000000000:1000000000:10000000000),(20000000000:10000000000:100000000000),(200000000000:100000000000:1000000000000)];
%variFewTslAsym=variAll(1:9:end);
%variAllCoV=0.015:0.005:0.50; 
%variFewCoV=[0.015,(0.05:0.05:0.50)];
%variAllXavg=100:100:9000; 
%variFewXavg=[100,(1000:1000:9000)];
%variAllk=[(1:0.005:1.09),(1.1:0.05:1.2),(1.25:0.05:1.9),(2:0.1:9),(10:1:100)]; 
%variFewk=[1,1.02,1.05,1.1,1.2,1.3,1.4,1.5,1.7,2,5,10,20,50,100]; 

variAll=100:100:9000; 
variFew=[100,(1000:1000:9000)];

Xjumpfac=100;
Xmaxfac=2;

numAll=length(variAll);
numFew=length(variFew);

statsAll=zeros(2*n+3,numAll+1);
statsAll(2:n+2,1)=(0:n)';
statsAll(n+3:2*n+3,1)=(0:n)';

statsFew=zeros(2*n+3,numFew+1);
statsFew(2:n+2,1)=(0:n)';
statsFew(n+3:2*n+3,1)=(0:n)';

iFew=1;
lnSuoFew0=zeros(1,numFew);
lnSuoFew1=zeros(1,numFew);
lnSuoFew2=zeros(1,numFew);
lnSuoFew5=zeros(1,numFew);
lnSuoFew10=zeros(1,numFew);
lnSuoFew15=zeros(1,numFew);
lnSuoFew20=zeros(1,numFew);
XFew=zeros(1,numFew);

for ivari=1:numAll
    clearvars -global -except ...
        Xavgo CoVo lnSuo X nX...
        ivari variAll numAll statsAll ...
        iFew variFew XFew nXFew statsFew Xjumpfac Xmaxfac...
        lnSuoFew0 lnSuoFew1 lnSuoFew2 lnSuoFew5 lnSuoFew10 lnSuoFew15 lnSuoFew20
    varinow=variAll(ivari);
    nominalX(n,varinow);
    statsAll(:,ivari+1)=[varinow, 1/1000*Xavgo, 100*CoVo]';
    if varinow==variFew(iFew)
        iXFew=(1:Xjumpfac:nX/Xmaxfac);
        XFewnow=X(iXFew);
        nXFew=length(XFewnow);
        statsFew(:,iFew+1)=statsAll(:,ivari+1);
        lnSuoFew0(1:nXFew+1,iFew)=[varinow;lnSuo(iXFew,1)];
        lnSuoFew1(1:nXFew+1,iFew)=[varinow;lnSuo(iXFew,2)];
        lnSuoFew2(1:nXFew+1,iFew)=[varinow;lnSuo(iXFew,3)];
        lnSuoFew5(1:nXFew+1,iFew)=[varinow;lnSuo(iXFew,6)];
        lnSuoFew10(1:nXFew+1,iFew)=[varinow;lnSuo(iXFew,11)];
        lnSuoFew15(1:nXFew+1,iFew)=[varinow;lnSuo(iXFew,16)];
       % lnSuoFew20(1:nXFew+1,iFew)=[varinow;lnSuo(iXFew,21)];
        XFew(1:length(XFewnow)+1,iFew)=[varinow;1/1000*XFewnow];
        iFew=iFew+1;
    end
    
end

statsAll=statsAll';

end

