%-----------------------------------------------------------
function uu=ushu51(x,ap)
% ex 4.5.: single solition
%c=0.3; x0=0.5; k=sqrt(c/ap)/2;
%uu=3*c*(sech(k*(x-x0))).^2;
% ex 4.5: double solition
c1=0.3; c2=0.1;x1=0.4;x2=0.8;
k1=0.5*sqrt(c1/ap); k2=0.5*sqrt(c2/ap);
uu=3*c1*(sech(k1*(x-x1))).^2+3*c2*(sech(k2*(x-x2))).^2;
% ex 4.5: triple soliton
%tmp=sqrt(108*ap);
%uu=2/3*(sech((x-1)/tmp)).^2; %.^ array power
% below for ex4.1. of Shu, so ap serves as t
%uu=sin(x+ap);

