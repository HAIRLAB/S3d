function  f=singlesd(x,t,c,a)

 f=c*sqrt(c)/4*sech(sqrt(c)/2*(x-c*t-a))*sech(sqrt(c)/2*(x-c*t-a))*sech(sqrt(c)/2*(x-c*t-a))*...
     (exp(-sqrt(c)/2*(x-c*t-a))-exp(sqrt(c)/2*(x-c*t-a)));