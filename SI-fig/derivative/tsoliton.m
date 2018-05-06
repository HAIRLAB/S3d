function u=tsoliton(x, t, a_1,a_2,s_1,s_2)

f_1=exp(-a_1*(x-s_1)+a_1*a_1*a_1*t);
f_2=exp(-a_2*(x-s_2)+a_2*a_2*a_2*t);
l=(a_2-a_1)/(a_2+a_1);
u=12*(a_1*a_1*f_1+a_2*a_2*f_2+2*(a_1-a_2)^2*f_1*f_2+l*l*(a_2*a_2*f_1*f_1*f_2+a_1*a_1*f_1*f_2*f_2))/(1+f_1+f_2+l*f_1*f_2)^2;
