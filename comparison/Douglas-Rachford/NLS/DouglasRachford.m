function x = DouglasRachford(A,b,sigma,tau,mu,MaxIt,tol)

% ============================================================
% References:
%   (1) Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
%   TBA
%   (2) Gabriel Peyre
%   http://www.numerical-tours.com/matlab/sparsity_5_dictionary_learning_denoising/
%
% Inputs:
%   A = m*n matrix
%   b = m*1 vector
%   sigma = constrain parameter
%   tau = step size
%   mu = convergence parameter
%   MaxIt = maximum number of iterations allowed
%   tol = tolerance
%
% Output:
%   u = min_{x}     ||x||_1
%       subject to  ||Ax-b||_2 <= sigma
%
% Authors: Hayden Schaeffer, Giang Tran, Rachel Ward, Linan Zhang
% Date: May 9, 2018
% ============================================================

% initialization
T = 1;
N = size(A,2);
M = length(b);
x = zeros(N,1);
x1 = zeros(N,1); % intermediate approximation for x
u1 = zeros(M,1); % intermediate approximation for an auxiliary variable
error = tol+1;

pAlower = chol(eye(N) + A'*A,'lower');
mu1 = 1-mu;

while T<=MaxIt && error>=tol
    
    % first step
    [p,q] = rProxF(x1,u1);
    [p,q] = rProxG(p,q);
    x1 = mu1*x1 + mu*p;
    u1 = mu1*u1 + mu*q;
    
    % second step
    [xnew,~] = ProxF(x1,u1);
    
    % update
    error = norm(x-xnew);
    x = xnew;
    T = T+1;
    
end

    function [x1,u1] = ProxF(x,u)
        x1 = max(abs(x)-tau,0) .* sign(x);
        u1 = b + (u-b) * min(sigma/norm(u-b,2),1);
    end

    function [x1,u1] = rProxF(x,u)
        [x1,u1] = ProxF(x,u);
        x1 = 2*x1-x;
        u1 = 2*u1-u;
    end

    function [x1,u1] = ProxG(x,u)
        y = x + A'*u;
        x1 = pAlower'\(pAlower\y);
        u1 = A*x1;
    end

    function [x1,u1] = rProxG(x,u)
        [x1,u1] = ProxG(x,u);
        x1 = 2*x1-x;
        u1 = 2*u1-u;
    end

end