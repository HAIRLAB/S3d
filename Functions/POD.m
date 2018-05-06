function [phi, lam, nbasis] = POD(X,cenergy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the POD basis vectors without
% subtracting the mean of the ensemble. Actually the ensemble
% mean is nothing more the first POD basis vector if we perform
% POD for the ensemble with nonzero mean.
%
% Function inputs:
% X       :  (n x nsnap) snapshot matrix
%            n = number of states in full space
%            nsnap = number of snapshots or  time partition
% cenergy :  how much energy of the ensemble you want 
%            to capture, i.e. cenergy = 99.9 (percent)
% 
% Function outputs:
% phi     :   (n x nsnap) matrix containing POD basis vectors
% lam     :   (nsnap x 1) vector containing POD eigenvalues
% nbasis  :   number of POD basis vectors you should use
%             to capture cenergy (percent) of the ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See also POD

%   METHOD OF SNAPSHOTS
% calculate the empirical correlation matrix C
C = X*X';
% Calculate the POD basis
% Note that the correlation matrix C is symmetric, so SVD and EIG
% will give the same evectorC and evalueC but they are already in
% descending order and hence we don't need to rearrange evectorC and evalueC
% if SVD is used
[leftmatrixC,singvalueC] = svd(C);
phi = leftmatrixC;
% return the POD eigenvalues in a vector
evalueC=singvalueC;
lam = diag(evalueC);
%%%       Find the number of POD basis vectors capturing cenergy (percent) of energy   %%%%
% total energy
tenergy = sum(lam);
energy = 0.;

i = 1;
while (((energy/tenergy)*100) < cenergy)
    energy = energy + lam(i);
    i = i+1;
end
nbasis = i;

% plot eigenvalues corresponding to nbasis vectors
%  plot(lam(1:nbasis)/tenergy/100,'*')
%  xlabel('Number of POD eigenvalues')
%  ylabel('Energy captured')