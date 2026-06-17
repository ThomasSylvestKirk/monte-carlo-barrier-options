function [Est,SE,HW]=BarrieroptionABAST(S0,D,K,sigma,T,n,repetetions,r,Z)
%This simulate the underlying asset using a simple simulation where the
%barriers are shifted.
h=T/n;
logS=nan(n+1,repetetions);
logS(1,:)=log(S0);
un=0.5826+0.1245*exp(-2.7*((log(S0/D))/(sigma*sqrt(h)))^1.2);
%Constant volatility
for j=1:repetetions
    hit=0;
for i=1:n
    logS(i+1,j)=logS(i,j)+(r-0.5*sigma^2)*h+sqrt(h)*sigma*Z(i,j);
    if logS(i+1,j)<log(D*exp(un*sigma*sqrt(h)));
        Y(j)=0;
    hit=1;
    break
    end
end
if hit==0
   Y(j)=exp(-r*T)*max(0,exp(logS(end,j))-K);
end
end
Est=mean(Y);
SE=sqrt(var(Y));
HW=SE*icdf('normal',[0.975],0,1);
end