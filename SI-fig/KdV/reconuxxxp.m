%---------------------------------------------------------
% reconstruction the 3rd derivative for periodic function
%---------------------------------------------------------
function ux=reconuxxxp(u,N,h)
IEX=0; % 4 means explicte 4th-order reconstruction
if IEX==4
a=2; b=-1;
for i=1:N
if i==1
tmp = b*(u(i+3)-3*u(i+1)+3*u(N-1)-u(N-3))/4 ...
+ a*(u(i+2)-2*u(i+1)+2*u(N-1)-u(N-2));
elseif i==2
tmp = b*(u(i+3)-3*u(i+1)+3*u(i-1)-u(N-2))/4 ...
+ a*(u(i+2)-2*u(i+1)+2*u(i-1)-u(N-1));
elseif i==3
tmp = b*(u(i+3)-3*u(i+1)+3*u(i-1)-u(N-1))/4 ...
+ a*(u(i+2)-2*u(i+1)+2*u(i-1)-u(i-2));
elseif i== (N-2)
tmp = b*(u(2)-3*u(i+1)+3*u(i-1)-u(i-3))/4 ...
+ a*(u(i+2)-2*u(i+1)+2*u(i-1)-u(i-2));
elseif i== (N-1)
    tmp = b*(u(3)-3*u(i+1)+3*u(i-1)-u(i-3))/4 ...
+ a*(u(2)-2*u(i+1)+2*u(i-1)-u(i-2));
elseif i==N
tmp = b*(u(4)-3*u(2)+3*u(i-1)-u(i-3))/4 ...
+ a*(u(3)-2*u(2)+2*u(i-1)-u(i-2));
else
tmp = b*(u(i+3)-3*u(i+1)+3*u(i-1)-u(i-3))/4 ...
+ a*(u(i+2)-2*u(i+1)+2*u(i-1)-u(i-2));
end
ux(i)=tmp/(2*h^3);
end
return
end
% 6 means explict 6th-order reconstruction
if IEX==6
a=-488/240; b=338/240; c=-72/240; d=7/240;
for i=1:N
if i==1
tmp = a*(u(i+1)-u(N-1))+b*(u(i+2)-u(N-2)) ...
+ c*(u(i+3)-u(N-3))+d*(u(i+4)-u(N-4));
elseif i==2
tmp = a*(u(i+1)-u(i-1))+b*(u(i+2)-u(N-1)) ...
+ c*(u(i+3)-u(N-2))+d*(u(i+4)-u(N-3));
elseif i==3
tmp = a*(u(i+1)-u(i-1))+b*(u(i+2)-u(i-2)) ...
+ c*(u(i+3)-u(N-1))+d*(u(i+4)-u(N-2));
elseif i==4
tmp = a*(u(i+1)-u(i-1))+b*(u(i+2)-u(i-2)) ...
+ c*(u(i+3)-u(i-3))+d*(u(i+4)-u(N-1));
elseif i==(N-3)
tmp = a*(u(i+1)-u(i-1))+b*(u(i+2)-u(i-2)) ...
+ c*(u(i+3)-u(i-3))+d*(u(2)-u(i-4));
elseif i== (N-2)
tmp = a*(u(i+1)-u(i-1))+b*(u(i+2)-u(i-2)) ...
+ c*(u(2)-u(i-3))+d*(u(3)-u(i-4));
elseif i== (N-1)
tmp = a*(u(i+1)-u(i-1))+b*(u(2)-u(i-2)) ...
+ c*(u(3)-u(i-3))+d*(u(4)-u(i-4));
elseif i==N
tmp = a*(u(2)-u(i-1))+b*(u(3)-u(i-2)) ...
+ c*(u(4)-u(i-3))+d*(u(5)-u(i-4));
else
tmp = a*(u(i+1)-u(i-1))+b*(u(i+2)-u(i-2)) ...
+ c*(u(i+3)-u(i-3))+d*(u(i+4)-u(i-4));
end
ux(i)=tmp/(h^3);
end
return
end
% below are implicit reconstruction
ISC=6; % order of the scheme
% for 4th-order compact scheme: alfa cannot be 1/2
if ISC==4
alfa=15/32; aa=2; bb=2*alfa-1;
% for 6th-order compact scheme
elseif ISC==6
alfa=7/16; aa=2; bb=-1.0/8;
end
amat=zeros(N,N);
B=zeros(N,1);
for i=1:N
if i==1
amat(1,1)=1; amat(1,2)=alfa; amat(1,N-1)=alfa;
B(i)=bb*(u(i+3)-3*u(i+1)+3*u(N-1)-u(N-3))/4 ...
+aa*(u(i+2)-2*u(i+1)+2*u(N-1)-u(N-2));
elseif i==2
amat(2,1)=alfa; amat(2,2)=1; amat(2,3)=alfa;
B(i)=bb*(u(i+3)-3*u(i+1)+3*u(i-1)-u(N-2))/4 ...
+aa*(u(i+2)-2*u(i+1)+2*u(i-1)-u(N-1));
elseif i==3
amat(i,i-1)=alfa; amat(i,i)=1; amat(i,i+1)=alfa;
B(i)=bb*(u(i+3)-3*u(i+1)+3*u(i-1)-u(N-1))/4 ...
+aa*(u(i+2)-2*u(i+1)+2*u(i-1)-u(i-2));
elseif i==N-2
amat(i,i-1)=alfa; amat(i,i)=1; amat(i,i+1)=alfa;
B(i)=bb*(u(2)-3*u(i+1)+3*u(i-1)-u(i-3))/4 ...
+aa*(u(i+2)-2*u(i+1)+2*u(i-1)-u(i-2));
elseif i==N-1
amat(i,i-1)=alfa; amat(i,i)=1; amat(i,i+1)=alfa;
B(i)=bb*(u(3)-3*u(i+1)+3*u(i-1)-u(i-3))/4 ...
+aa*(u(2)-2*u(i+1)+2*u(i-1)-u(i-2));
elseif i==N
amat(N,1)=-1; amat(N,N)=1;
B(i)=0;
else % i>=4 & i <=N-3
    amat(i,i-1)=alfa; amat(i,i)=1; amat(i,i+1)=alfa;
B(i)=bb*(u(i+3)-3*u(i+1)+3*u(i-1)-u(i-3))/4 ...
+aa*(u(i+2)-2*u(i+1)+2*u(i-1)-u(i-2));
end
B(i)=B(i)/(2*h^3);
end
% call trisys.m
ux=(amat\B)';