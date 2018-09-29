function x = DouglasRachfordxt(A,b,ga,lammda,mu,MaxIt,tol)

% ============================================================
%
% Inputs:
%   A = m*n matrix
%   b = m*1 vector
%   ga = constrain parameter   1/2
%   lammda = balance parameter   500
%   mu = convergence parameter  0.5
%   MaxIt = maximum number of iterations allowed  5000
%   tol = tolerance
% 
% Output:
%   u = min_{x}     ||x||_1
%       subject to  ||Ax-b||_2 <= sigma
%
% Authors: Hayden Schaeffer
% ============================================================

% initialization
T = 1;
N = size(A,2);
M = length(b);
x = zeros(N,1);
x1 = zeros(N,1); % intermediate approximation for x
error = tol+1;

pAlower = eye(N) + ga.*A'*A;
mu1 = 1-mu;

while T<=MaxIt && error>=tol
    
    % first step
    [p] = rProxF(x1); 
    [p] = rProxG(p);
    x1 = mu1*x1 + mu*p; 
    
    % second step
    [xnew] = ProxF(x1);
    
    % update
    error = norm(x-xnew);
    x = xnew;
    T = T+1;
    
end

    function [x1] = ProxF(x)   %%F--H1
        x1 = max(abs(x)-ga.*lammda,0) .* sign(x);
    end

    function [x1] = rProxF(x)
        [x1] = ProxF(x);
        x1 = 2*x1-x;
    end

    function [x1] = ProxG(x)  %G--H2
        y = x + ga.*A'*b;
        x1 = pinv(pAlower)*y;
    end

    function [x1] = rProxG(x)
        [x1] = ProxG(x);
        x1 = 2*x1-x;
    end

end