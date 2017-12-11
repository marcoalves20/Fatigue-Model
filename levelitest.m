function [lnSurI,lnSlrI,lnSkrI,lnSueI,lnSkeI,lnSueI1] =...
    leveli(lnSueO,lnSkeO,lnSueO1,AI,CI,leO,lef,lcontrol,deltaA,lintact,nXcritical,lambda,Tf,tau0)

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
%--------------------------------------------------------------------------
%Define auxiliar parameters, to seperate the different versions of the
%scalling law

for i=1:nX
    if lef(i)>=Lout/2
        nAux1 = i;
        break
    end
end

for i=1:nX
    if deltaA(i)>=Lout
        nXcritical = i;
        break
    end
end
%-------------------------------------------------------------------------
if nAux1 ~= nXcritical
    %Original Scalling law
    if k==1
      lnSucI(1:nAux1,1)=2*lnSueO(1:nAux1,1)+log(2-exp(2*lnSueO(1:nAux1,1)));
    else
      if Ck > 3
        lnSucI(1:nAux1,1)=4*lnSueO(1:nAux1,1)+log(1+2*exp(lnSkeO(1:nAux1,1)-3*lnSueO(1:nAux1,1))-2*exp(lnSkeO(1:nAux1,1)-lnSueO(1:nAux1,1)));
   	    lnSucI(lnSucI==Inf)=lnSueO(lnSucI==Inf)+lnSkeO(lnSucI==Inf)+log(2+exp(3*lnSueO(lnSucI==Inf)-lnSkeO(lnSucI==Inf))-2*exp(2*lnSueO(lnSucI==Inf)));      
      else
        lnSucI(1:nAux1,1)=lnSueO(1:nAux1,1)+lnSkeO(1:nAux1,1)+log(2+exp(3*lnSueO(1:nAux1,1)-lnSkeO(1:nAux1,1))-2*exp(2*lnSueO(1:nAux1,1))); 
        lnSucI(lnSucI==Inf)=4*lnSueO(lnSucI==Inf)+log(1+2*exp(lnSkeO(lnSucI==Inf)-3*lnSueO(lnSucI==Inf))-2*exp(lnSkeO(lnSucI==Inf)-lnSueO(lnSucI==Inf)));
      end
    end

    %Scalling law changes because lc ~= 2*le
    lnSucI(nAux1+1:nXcritical)=2*lnSueO(nAux1+1:nXcritical)+2*lnSueO1(nAux1+1:nXcritical)+log(1+2*exp(lnSkeO(nAux1+1:nXcritical)-2*lnSueO(nAux1+1:nXcritical)-lnSueO1(nAux1+1:nXcritical))-2*exp(lnSkeO(nAux1+1:nXcritical)-lnSueO1(nAux1+1:nXcritical)));
    %lnSucI(lnSucI==Inf)=lnSkeO(lnSucI==Inf)+log(exp(4*lnSueO(lnSucI==Inf)-lnSkeO(lnSucI==Inf))+2*exp(lnSueO1(lnSucI==Inf))-2*exp(2*lnSueO(lnSucI==Inf)+lnSueO1(lnSucI==Inf)));

    % Modified Scalling law
    lnSucI(nXcritical+1:nX,1)=2*lnSueO(nXcritical+1:nX,1)+log(1+2*exp(lnSkeO(nXcritical+1:nX,1)-2*lnSueO(nXcritical+1:nX,1))-2*exp(lnSkeO(nXcritical+1:nX,1)-lnSueO(nXcritical+1:nX,1)));
    lnSucI(lnSucI==Inf)=lnSkeO(lnSucI==Inf)+log(2+exp(2*lnSueO(lnSucI==Inf)-lnSkeO(lnSucI==Inf))-2*exp(lnSueO(lnSucI==Inf)));
else 
    if k==1
      lnSucI(1:nAux1,1)=2*lnSueO(1:nAux1,1)+log(2-exp(2*lnSueO(1:nAux1,1)));
    else
      if Ck > 3
         lnSucI(1:nAux1,1)=4*lnSueO(1:nAux1,1)+log(1+2*exp(lnSkeO(1:nAux1,1)-3*lnSueO(1:nAux1,1))-2*exp(lnSkeO(1:nAux1,1)-lnSueO(1:nAux1,1)));
   	    lnSucI(lnSucI==Inf)=lnSueO(lnSucI==Inf)+lnSkeO(lnSucI==Inf)+log(2+exp(3*lnSueO(lnSucI==Inf)-lnSkeO(lnSucI==Inf))-2*exp(2*lnSueO(lnSucI==Inf)));      
      else
            lnSucI(1:nAux1,1)=lnSueO(1:nAux1,1)+lnSkeO(1:nAux1,1)+log(2+exp(3*lnSueO(1:nAux1,1)-lnSkeO(1:nAux1,1))-2*exp(2*lnSueO(1:nAux1,1))); 
            lnSucI(lnSucI==Inf)=4*lnSueO(lnSucI==Inf)+log(1+2*exp(lnSkeO(lnSucI==Inf)-3*lnSueO(lnSucI==Inf))-2*exp(lnSkeO(lnSucI==Inf)-lnSueO(lnSucI==Inf)));
      end
    end
    
     % WLT is applied for stresses nXcrit:nX
      % lnSucI(nXcritical+1:nX,1)=4*lnSueO(nXcritical+1:nX,1);

    % Modified Scalling law
    lnSucI(nXcritical+1:nX,1)=2*lnSueO(nXcritical+1:nX,1)+log(1+2*exp(lnSkeO(nXcritical+1:nX,1)-2*lnSueO(nXcritical+1:nX,1))-2*exp(lnSkeO(nXcritical+1:nX,1)-lnSueO(nXcritical+1:nX,1)));
    lnSucI(lnSucI==Inf)=lnSkeO(lnSucI==Inf)+log(2+exp(2*lnSueO(lnSucI==Inf)-lnSkeO(lnSucI==Inf))-2*exp(lnSueO(lnSucI==Inf)));
end
%--------------------------------------------------------------------------

lnSucI(1)=0;
%Scalling back to the reference length-------------------------------------
lnSurI=zeros(nX,1);
% lnSurI(:,1)=1./(2*leO).*lnSucI; --> WHY?
lnSurI(:,1)=1./(lcontrol).*lnSucI;
lnSurI(1)=0;
%-------------------------------------------------------------------------

lnSkrI1=zeros(nX,1);
lnSkrI2=zeros(nX,1);
lnSkeI1=zeros(nX,1);
lnSkeI2=zeros(nX,1);
%Calculates SkrI1 up to X and then it is limited to nXcrit in SkrI1
 lnSlrI=cumtrapz(lnSurI(1:nX))*DX./X;
 lnSlrI(1)=0;
    if mod(k,1)==0
        %If k is integer, takes every other k-value:
        lnSkrI1(1:nK,1)=(k*lnSlrI(1:k:end)-lnSlrI(1:nK))./(k-1);
    else
        %If k is not integer (interpolation required):   
        %Interpolates lnSlr for k*X (non-integer k):
        lnSlrKI(1:nK,1)=interp1(X,lnSlrI,XK,'linear');
        %uses lnSlr interpolated @ k*X (lnSlrK)
        lnSkrI1(1:nK,1)=(k*lnSlrKI-lnSlrI(1:nK))./(k-1);
    end
    lnSkrI1(nK+1:nX,1)=Ck*lnSurI(nK+1:nX);
    lnSkrI1(1)=0;
    %--------------------------------------------------------------------------  
    %Survival probability calculated for the stress plateau
    %Modification 21/03/17----------------------------------
    lnSkrI2(1:nK)=lnSurI(1:k:end);
    lnSkrI2(nK+1:end)=k^m*lnSurI(nK+1:end);
    lnSkrI2(1)=0;
    %-------------------------------------------------------
    lnSkrI=lnSkrI2+lnSkrI1;
    %%%For effective recovery length%%%%%%%%%%%%%%%%%%
    %Calculates survival functions    
    lnSkeI1=lintact.*lnSkrI1;
    lnSkeI2=deltaA.*lnSkrI2;
    lnSkeI=lnSkeI1+lnSkeI2;
    lnSueI=(lef).*lnSurI;
    lnSueI1=(lcontrol-lef).*lnSurI;

end

