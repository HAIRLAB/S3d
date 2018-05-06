function rhs=nls_rhs(t,ut,dummy,k)
u=ifft(ut);
rhs=-(3*i/10)*(k.^2).*ut+(10*i/3)*fft( (abs(u).^2).*u -u);