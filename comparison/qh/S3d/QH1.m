% 
% @author: Ye Yuan
% 
% utilizing S3d method to analyze the QH equation:            
%          u{t} = 0.5iu{xx} - i(x^2/2)u
%

addpath ../../../Functions

load('data1xt500.mat');

theta1 = [real(theta1) -imag(theta1);imag(theta1) real(theta1) ];
y1 = double([real(y1);imag(y1)]);

%normalization
T1 = zeros(size(theta1,1),size(theta1,2));
Mreg1 = zeros(size(theta1,2),1);
for i = 1:size(theta1,2)
    Mreg1(i,1) = 1/max(abs(theta1(:,i)));
    T1(:,i) = Mreg1(i)*theta1(:,i);    
end

MAXITER = 7;    %maximum number of inner iterations
lambda1 = 0.25;  %the regularization parameter 
w1 = tac_reconstruction(y1, T1, lambda1,MAXITER);
wt1 = w1(:,1:end).*Mreg1; %identified parameters corresponding to the basis function in Theta

w1 = w1(:,end).*Mreg1; %identified parameters corresponding to the basis function in Theta

%print result and error

threshold = 1e-6; %filter those cofficients in w being less than the threshold
y_name = 'u{t}';

for j= 1:size(wt1,2)
for i = 1:size(wt1,1)
    if abs(wt1(i,j))<threshold
        wt1(i,j) = 0;
    end
end
end

wt = wt1(1:40,:) + sqrt(-1)*wt1(41:80,:);

