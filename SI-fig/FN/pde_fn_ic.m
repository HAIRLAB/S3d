function u0 = pde_fn_ic(x,A,B,C)
u0 = [1*exp(-((x)/1).^2)
0.2*exp(-((x+2)/1).^2)];

% u0=[sin(pi*x*0.5)
%     sin(pi*x*0.5)];

%  u0=[5*exp(-0.5*x.^2)
%      zeros(size(x))];
