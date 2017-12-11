clear all
clc
% x=[1 2 3 4 5];
% y=[-3 4 1 5 3];
% z=[2 3 5 8 2];
% plot(x,y,'Color',[0 0 1],'Marker','+')
% hold on
% plot(x,z,'Color',[1 0 1],'Marker','+');

x = linspace (0 ,2*pi ,100) ;
y = sin( x ) ;
plot (x ,y ,'Color', 1/255*[0 205 0]) ;
axis ([0 2* pi -1 1]) ;
title ('sin (\ theta )',' FontSize ' ,18) ;
xlabel ('\ theta ',' FontSize ' ,16) ;
