%---------------------------------------------------------
% implement 10th order filter for periodic function
%---------------------------------------------------------
function uf=filterLI(u,N)
ID=5; %10th-order filter
af=0.1;
a(1)=(193+126*af)/256;
a(2)=(105+302*af)/256;
a(3)=15*(-1+2*af)/64;
a(4)=45*(1-2*af)/512;
a(5)=5*(-1+2*af)/256;
a(6)=(1-2*af)/512;
amat=zeros(N,N);
B=zeros(N,1);
for i=1:N
if i==1
amat(1,1)=1; amat(1,2)=af; amat(1,N-1)=af;
elseif i==N
amat(N,1)=-1; amat(N,N)=1;
else % special for the 1st and last rows
amat(i,i-1)=af; amat(i,i)=1; amat(i,i+1)=af;
end
B(i)=0;
for k=0:ID
% try to keep i+k and i-k into the range of [1,N]
k1=i+k;
if k1 > N
k1=mod(k1,N-1);
end
k2=i-k;
if k2 <=0
k2=k2+(N-1);
end
%disp(k1), disp(k2),
B(i)=B(i)+a(k+1)*(u(k1)+u(k2));
end
%disp(¡¯----------¡¯);
B(i)=B(i)/2;
end
B(N)=0; % the last row is special
% call linear system solver
uf=(amat\B)';