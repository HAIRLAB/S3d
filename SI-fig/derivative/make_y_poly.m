function [ ut] = make_y_poly( u , t,t_deg,parameter )

x_cut = parameter.x_cut ;
t_cut = parameter.t_cut ;
order = parameter.order ;
x_num_to_fit = parameter.x_num_to_fit ;

[r,v] = size(u);
ut = zeros(r-2*x_cut,v-2*t_cut);
num = (r-2*x_cut)*(v-2*t_cut);
all_derivative.derivative = reshape(u(1+x_cut:end-x_cut,1+t_cut:end-t_cut),1,num)';


Fit = zeros(size(ut,1)*size(ut,2),order+1);
for i = 1: size(ut,2)
    for k = 1:size(ut,1)
        Fit((i-1)*size(ut,1)+k,:) = polyfit(x(k : k+2*x_num_to_fit)',u(k : k+2*x_num_to_fit,i+t_cut),order);
    end
end




for i = 1: size(ut,1)
    for k = 1:size(ut,2)
        tmp = polyfit(t(k:k+2*t_num_to_fit)',u(i+x_cut,k:k+2*t_num_to_fit),order);
        for m = 1:t_deg
            tmp = polyder(tmp );
        end
        ut(i,k) = polyval(tmp,t(k+t_num_to_fit));
    end
end
ut = reshape(ut,1,num)';

end

