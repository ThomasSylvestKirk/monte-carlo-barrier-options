function p=CSP(X,n,D,sigma,monitoringinterval)
%Calculate the conditional survival probability. It uses value of the
%underlying asset on every monitoring points.
p=1;
if X>D
for i=1:n 
        value=-2*((X(i)-D)*(X(i+1)-D))/(sigma^2*monitoringinterval);
   if value>-36.044
        p=p*(1-exp(value));
    end
end
end
end