function [c,f,s] = pde_fn_pde(x,t,u,DuDx,A,B,C)
c = [1; 1]; % c diagonal terms
f = [1; 0] .* DuDx; % diffusion term
rhsV = u(1)*(A-u(1))*(u(1)-1) - u(2);
rhsW = B*u(1)-C*u(2);
s = [rhsV; rhsW];