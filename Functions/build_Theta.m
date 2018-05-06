function [theta,theta_name] = build_Theta(data,derdata,dataname,dername,polyorder)

% @author: Ye Yuan
%This function builds the candidate dictionaries that represent 
%polynoimials up to degree polyorder of all input variables
%
%Inputs:
%      data      : sample data of states
%      derdata   : data of the derivative terms
%      dataname  : names of variables in data
%      dername   : names of derivatives in derdata
%      polyorder : maximum power of polynomial
%
%Outputs:
%      theta     : dictionary matrix 
%      theta_name: names of variables in theta 
%

theta(:,1) = ones(size(data,1),1);
theta_name{1} = '1';

theta_name = [theta_name,dataname];
theta = [theta,data];

%generate all polynomials up to degree polyorder of variables in data 
if polyorder > 1
    for i = 2:polyorder
        tmp = ones(1,i);
        comb = [];        
        [comb] = comb_with_rep(comb,tmp,size(data,2),i); %generate all possible combinations
        for j = 1:size(comb,1)
            nametmp = dataname(comb(j,1));
            datatmp = data(:,comb(j,1));
            for k = 2:size(comb,2) 
                nametmp = strcat(nametmp,'*',dataname(comb(j,k)));
                datatmp = datatmp.*data(:,comb(j,k));
            end
            theta = [theta,datatmp];
            theta_name = [theta_name,nametmp];
        end
    end
end

%add the derivatives into theta 
len = length(theta_name);
theta = [theta,derdata];
theta_name = [theta_name,dername];

%add polynomials times derivatives into theta
for i = 2:len
    for j = 1:length(dername)
        datatmp = theta(:,i).*derdata(:,j);
        nametmp = [theta_name{i} '*' dername{j}];
        theta = [theta,datatmp];
        theta_name = [theta_name,nametmp];
    end
end

end



function [comb] = comb_with_rep(comb,tmp,N,K)

% @author: Ye Yuan
%This recursive function generates combinations with replacement.
%
%Inputs:
%      comb: combinations that have been created
%      tmp : next combination to be included in comb
%      N   : number of all possible elements
%      K   : number of elements to be selected 
%
%Output:
%      comb: updated combinations


if tmp(1)>N
    return;
end

for i = tmp(K):N
    comb = [comb;tmp];
    tmp(K) = tmp(K) + 1;
end
tmp(K) = tmp(K) - 1;
for i = K:-1:1
    if tmp(i)~=N
        break;
    end
end
tmp(i) = tmp(i)+1;
for j = i+1:K
    tmp(j)=tmp(i);
end

[comb] = comb_with_rep(comb,tmp,N,K);

end
