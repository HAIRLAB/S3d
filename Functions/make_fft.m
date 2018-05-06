function [ut] = make_fft(u,order,L)

% @author: Ye Yuan
%This function computes the derivative of given data by using Fourier 
%spectral method to approximate derivative
%
%Inputs:
%      u    : one-dimensional data
%      order: order of derivative
%      L    : length of a period
%
%Output:
%      ut   : order-th derivative of u
%

n = size(u,1);

switch mod(n,2)
    case 1
        k=(2*pi/L)*[0:(n/2-1) 0  (-n/2):-1]';
        
    case 0
        k=(2*pi/L)*[0:n/2-1 0 -n/2+1:-1]';
end
g = 1i*k;

tmp = fft(u); 
ut = (g.^order).*tmp;
ut = real(ifft(ut));

end

