function [EST,SE,HW,Y]=DOCBB(S0,D,K,sigma,T,n,repetetions,r,Z)
%ONLY USES CONTINUOUS VALUE. This function estimate the price of the DOC and each of the expected payoffs which can be used to estiamte the central difference
S=nan(n+1,repetetions);
logs=nan(n+1,repetetions);

PDISK=nan(repetetions,1);
S(1,:)=S0;
logs(1,:)=log(S0);
Y=nan(repetetions,1);
h=T/n;
LD=log(D);
for j=1:repetetions
    p=1; hit=0; PDISK=1;
    for i=1:n
        logs(i+1,j)=logs(i,j)+(h*(r-sigma^2/2)+sigma*Z(i,j)*sqrt(h));
    end
    PDISK=1-P(logs(1,j),logs(end,j),h,LD,sigma,n);
    Y(j)=exp(-r*T)*PDISK*max(exp(logs(end,j))-K,0);
end
EST=mean(Y);
SE=std(Y);
HW=SE*icdf('normal',[0.975],0,1);
end