function [Estimate,SE,HW,Y]=BarrieroptionEulerNAIVwithout(S0,D,K,sigma,T,n,repetetions,r,Z)
%This simulate the underlying asset with a euler discreterization
h=T/n;
logS=nan(n+1,repetetions);
logS(1,:)=log(S0);
logD=log(D);
for j=1:repetetions
    hit=0;
for i=1:n
    logS(i+1,j)=logS(i,j)+(r-0.5*sigma^2)*h+sqrt(h)*sigma*Z(i,j); % Simulate the asset by Euler
    if logS(i+1,j)<logD
    Y(j)=0;
    hit=1;
    break
    end
end
if hit==0
   Y(j)=exp(-r*T)*max(0,exp(logS(end,j))-K);
end
end
Estimate=mean(Y);
SE=sqrt(var(Y));
HW=SE*icdf('normal',[0.975],0,1);
end