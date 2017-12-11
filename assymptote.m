function [  ] = assymptote( n )
%Assymtotic analysis
%1)WLT test

global lnSuo
global nX nK
global testWLT testXavgoWLT X DX

testWLT(:,1:n)=lnSuo(:,1:n);
for i=2:n
    testWLT(:,i)=2*testWLT(:,i-1);
end

for i=1:n
    testFuoWLT(:,i)=1-exp(testWLT(:,i));
    intFuoDX=trapz(testFuoWLT(:,i))*DX;
    testXavgoWLT(i)=X(end)-intFuoDX;
end

%2) Modified scalling law - Sa+b = Sa(sig)*Sb(sig)+2*(1-Sa(sig))*Sb(2sig)

% testWLT(1:nK,1:n)=lnSuo(1:nK,1:n);
% for i=2:n
% %     testWLT(1:nK,i)=2*(1-exp(testWLT(1:nK,i-1))).*exp(lnSuo(1:2:nX,i-1));
%     testWLT(1:nK,i)=exp(2*lnSuo(1:2:nX,i-1) - 2* testWLT(1:nK,i-1).*(lnSuo(1:2:nX,i-1)));
% end
% 
% for i=1:n
%     testFuoWLT(:,i)=1-testWLT(:,i);
%     intFuoDX=trapz(testFuoWLT(:,i))*DX;
%     testXavgoWLT(i)=X(end)-intFuoDX;
% end

end

