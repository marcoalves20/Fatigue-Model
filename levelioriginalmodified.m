function [lnSurI,lnSlrI,lnSkrI,lnSueI,lnSkeI] =...
    leveli(lnSueO,lnSkeO,AI,CI,leO,lef,deltaA,lintact,nXcritical,lambda,Tf,tau0)
%runs bundle strength model for generic level i
%called by bundleX

%%%Declaring variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
global DX
global Tsl
global Lref Lout
global k m
%Internal
global Ck 
global nX nK XK 
%Output
global X Xref

%%%Scaling law: calculates survival function in the control length%%%%%%%%%
%Analytical
if nX < nXcritical
        if k==1
        lnSucI=2*lnSueO+log(2-exp(2*lnSueO));
    else
        if Ck>3
            lnSucI=4*lnSueO+log(1+2*exp(lnSkeO-3*lnSueO)-2*exp(lnSkeO-lnSueO));
            lnSucI(lnSucI==Inf)=lnSueO(lnSucI==Inf)+lnSkeO(lnSucI==Inf)+log(2+exp(3*lnSueO(lnSucI==Inf)-lnSkeO(lnSucI==Inf))-2*exp(2*lnSueO(lnSucI==Inf)));
        else
           lnSucI=lnSueO+lnSkeO+log(2+exp(3*lnSueO-lnSkeO)-2*exp(2*lnSueO)); 
           lnSucI(lnSucI==Inf)=4*lnSueO(lnSucI==Inf)+log(1+2*exp(lnSkeO(lnSucI==Inf)-3*lnSueO(lnSucI==Inf))-2*exp(lnSkeO(lnSucI==Inf)-lnSueO(lnSucI==Inf)));
        end
    end
else
        if k==1
        lnSucI(1:nXcritical,1)=2*lnSueO(1:nXcritical,1)+log(2-exp(2*lnSueO(1:nXcritical,1)));
        else
        if Ck>3
            lnSucI(1:nXcritical,1)=4*lnSueO(1:nXcritical,1)+log(1+2*exp(lnSkeO(1:nXcritical,1)-3*lnSueO(1:nXcritical,1))-2*exp(lnSkeO(1:nXcritical,1)-lnSueO(1:nXcritical,1)));
            lnSucI(lnSucI==Inf)=lnSueO(lnSucI==Inf)+lnSkeO(lnSucI==Inf)+log(2+exp(3*lnSueO(lnSucI==Inf)-lnSkeO(lnSucI==Inf))-2*exp(2*lnSueO(lnSucI==Inf)));       
        else
           lnSucI(1:nXcritical,1)=lnSueO(1:nXcritical,1)+lnSkeO(1:nXcritical,1)+log(2+exp(3*lnSueO(1:nXcritical,1)-lnSkeO(1:nXcritical,1))-2*exp(2*lnSueO(1:nXcritical,1))); 
           lnSucI(lnSucI==Inf)=4*lnSueO(lnSucI==Inf)+log(1+2*exp(lnSkeO(lnSucI==Inf)-3*lnSueO(lnSucI==Inf))-2*exp(lnSkeO(lnSucI==Inf)-lnSueO(lnSucI==Inf)));
        end
        end
    % WLT is applied for stresses nXcrit:nX
          lnSucI(nXcritical+1:nX,1)=4*lnSueO(nXcritical+1:nX,1);
end

% -----------------------------------------------------------------------------
%%%For the reference length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates uniform survival function
lnSurI(:,1)=1./(2*leO).*lnSucI;
lnSurI(1)=0;
if nX<nXcritical
%Calculates linear survival function
lnSlrI=cumtrapz(lnSurI)*DX./X;
lnSlrI(1)=0;
%Calculates concentration survival function
      if mod(k,1)==0
          %If k is integer, takes every other k-value:
          lnSkrI(1:nK,1)=(k*lnSlrI(1:k:end)-lnSlrI(1:nK))./(k-1);
      else
          %If k is not integer (interpolation required):   
          %Interpolates lnSlr for k*X (non-integer k):
          lnSlrKI(1:nK,1)=interp1(X,lnSlrI,XK,'linear');
          %uses lnSlr interpolated @ k*X (lnSlrK)
          lnSkrI(1:nK,1)=(k*lnSlrKI-lnSlrI(1:nK))./(k-1);
      end
      
      lnSkrI(nK+1:nX,1)=Ck*lnSurI(nK+1:nX,1);
      %lnSkrI(nK+1:nX,1)=NaN;
      %lnSkrI(nK+1:nX,1)=Ck*lnSurI(nK+1:nX)-1./leO(nK+1:nX,1)*log(2)/2*(Ck-log(k)/(k-1));
      lnSkrI(1)=0;
      lnSkeI=(lef).*lnSkrI;
      lnSueI=lef.*lnSurI;
else
%-------------------------------------------------------------------------
%Calculates SkrI1 up to X and then it is limited to nXcrit in SkrI1
 lnSlrI=cumtrapz(lnSurI(1:nX))*DX./X;
 lnSlrI(1)=0;
    if mod(k,1)==0
        %If k is integer, takes every other k-value:
        for i = 1:nK
            lnSkrI1(i,1)=1/(X(i))*(trapz(lnSurI(i:(2*i-1)))*DX);
        end
    else
        %If k is not integer (interpolation required):   
        %Interpolates lnSlr for k*X (non-integer k):
        lnSlrKI(1:nK,1)=interp1(X,lnSlrI,XK,'linear');
        %uses lnSlr interpolated @ k*X (lnSlrK)
        lnSkrI1(1:nK,1)=(k*lnSlrKI-lnSlrI(1:nK))./(k-1);
    end
    lnSkrI1(nK+1:nX,1)=Ck*lnSurI(nK+1:nX);
%     lnSkrI1(nK+1:nX,1)=0;
    lnSkrI1(1)=0;
    %--------------------------------------------------------------------------  
    %Survival probability calculated for the stress plateau
    lnSkrI2=zeros(nX,1);
    lnSkeI1=zeros(nX,1);
    lnSkeI2=zeros(nX,1);
    %Modification 21/03/17----------------------------------
    lnSkrI2(1:nK)=lnSurI(1:k:end);
    lnSkrI2(nK+1:end)=k^m*lnSurI(nK+1:end);
    %-------------------------------------------------------
    lnSkrI=lnSkrI2+lnSkrI1; 
    %%%For effective recovery length%%%%%%%%%%%%%%%%%%
    %Calculates survival functions    
    lnSkeI1=lintact.*lnSkrI1;
    lnSkeI2=deltaA.*lnSkrI2;
    lnSkeI=lnSkeI1+lnSkeI2;
    lnSueI=lef.*lnSurI;

% lnSkeI(nXcrit+1:end)=deltaA(nXcrit+1:end).*lnSkrI2(nXcrit+1:end);
end
end

