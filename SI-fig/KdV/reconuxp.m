% Reconstruction for periodic function
function ux=reconuxp(u,N,h)
IEX=0; % 1 for explicit 4th-order scheme
if IEX==1
alfa=0; aa=2/3*(alfa+2); bb=(4*alfa-1)/3;
for i=1:N
if i==1
tmp =bb*(u(i+2)-u(N-2))/2+aa*(u(i+1)-u(N-1));
elseif i==2
tmp =bb*(u(i+2)-u(N-1))/2+aa*(u(i+1)-u(i-1));
elseif i==(N-1)
tmp=bb*(u(2)-u(i-2))/2+aa*(u(i+1)-u(i-1));
elseif i==N
tmp=bb*(u(3)-u(i-2))/2+aa*(u(2)-u(i-1));
else
tmp=bb*(u(i+2)-u(i-2))/2+aa*(u(i+1)-u(i-1));
end
ux(i)=tmp/(2*h);
end
return
end
% below are implicit scheme
ISC=6; % order of the scheme
% for C4 scheme
if ISC==4
alfa=1.0/4.0; aa=3.0/2; bb=0;
% for C6 scheme
elseif ISC==6
alfa=1.0/3; aa=14.0/9; bb=1.0/9;
end
amat=zeros(N,N);
B=zeros(N,1);
for i=1:N
if i==1
    amat(1,1)=1; amat(1,2)=alfa; amat(1,N-1)=alfa;
B(i)=bb*(u(3)-u(N-2))/4+aa*(u(2)-u(N-1))/2;
elseif i==2
amat(2,1)=alfa; amat(2,2)=1; amat(2,3)=alfa;
B(i)=bb*(u(4)-u(N-1))/4+aa*(u(3)-u(1))/2;
elseif i==N-1
amat(i,i-1)=alfa; amat(i,i)=1; amat(i,i+1)=alfa;
B(i)=bb*(u(2)-u(N-3))/4+aa*(u(N)-u(N-2))/2;
elseif i==N
amat(N,1)=-1; amat(N,N)=1;
B(i)=0;
else % i>=3 & i <=N-2
amat(i,i-1)=alfa; amat(i,i)=1; amat(i,i+1)=alfa;
B(i)=bb*(u(i+2)-u(i-2))/4+aa*(u(i+1)-u(i-1))/2;
end
B(i)=B(i)/h;
end
% call trisys.m
ux=(amat\B)';