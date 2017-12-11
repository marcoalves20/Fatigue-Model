function [] = WLT(n)
%runs nominal model and calculates WLT approximation

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global DX X nX
%Output for 
global lnSuo
global FuoWLT XavgoWLT CoVoWLT lnSuoWLT

%%%Running nominal model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fatigueX(n,2);

%Preallocating matrices
lnSuoWLT=zeros(nX,n+1,n+1);

%%%Uses WLT to calculate strength distributions, from all levels %%%%%%%%%%
for i=1:n+1
    lnSuoWLT(:,1:i,i)=lnSuo(:,1:i);
    for j=i+1:n+1
        lnSuoWLT(:,j,i)=2*lnSuoWLT(:,j-1,i);
    end
    
    FuoWLT(:,:,i)=1-exp(lnSuoWLT(:,:,i)); 
    
    intFuoDX=trapz(FuoWLT(:,:,i))*DX;
    XavgoWLT(i,:)=X(end)-intFuoDX;

    intXFuoDX=zeros(1,n+1);
    for j=1:n+1
        intXFuoDX(1,j)=trapz(X.*FuoWLT(:,j,i))*DX;
    end

    CoVoWLT(i,:)=sqrt((X(end)-XavgoWLT(i,:)).^2-2.*intXFuoDX+2.*XavgoWLT(i,:).*intFuoDX)./XavgoWLT(i,:);
    
end
