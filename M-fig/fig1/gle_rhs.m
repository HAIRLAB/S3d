function rhs=gle_rhs(t,uvt,dummy,rbeta0,rbeta1,rbeta2,rbeta3,ralp,lbeta0,lbeta1,lbeta2,lbeta3,lalp,k,n)
u=ifft(uvt(1:n));
v=ifft(uvt(n+1:2*n));

% Reaction Terms
u3=abs(u).^2.*u; v3=abs(v).^2.*v; u1v2=(abs(v).^2).*u; u2v1=(abs(u).^2).*v;
utrhs=fft(rbeta1*u+rbeta2*u3+rbeta3*u1v2);
vtrhs=fft2(lbeta1*v+lbeta2*v3+lbeta3*u2v1);

rhs=[-1i*ralp*k.*uvt(1:n)-rbeta0*(k.^2).*uvt(1:n)+utrhs
     1i*lalp*k.*uvt(n+1:2*n)-lbeta0*(k.^2).*uvt(n+1:2*n)+vtrhs];
 
 
