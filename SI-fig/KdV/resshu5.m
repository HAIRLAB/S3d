% work on kdv equ of Shu¡¯s ex4.5
% the residual R=-0.5*(u^2)_x-ap*(u_xxx)
% i.e., u1=(u^2)_x, u2=u_xxx
function uu=resshu5(u1,u2,ap)
uu=-0.5*u1-ap*u2;