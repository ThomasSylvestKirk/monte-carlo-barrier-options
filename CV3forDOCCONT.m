function [EST,SE,HW,b] = CV3forDOCCONT(S0,D,K,sigma,T,n,repetetions,r,Z)
%This functions simulates the option and use regress to estimate b. Uses control
%variates to estimate the payoff.
h=T/n;
repetetionsforest=repetetions/10;
[Y,X]=DOC3CV(S0,D,K,sigma,T,n,repetetionsforest,r,Z(:,1:repetetionsforest));
DC=D*exp(-0.5826*sigma*sqrt(h));
E= [1;DOCbls(S0,T,K,r,sigma,DC);blsprice(S0,K,r,T,sigma);S0];
b=regress(Y',X');
clear Y X
repetetions=repetetions-repetetionsforest;
[Y,X]=DOC3CV(S0,D,K,sigma,T,n,repetetions,r,Z(:,repetetionsforest+1:end));
for i=1+1:length(b)
Y=Y-b(i)*(X(i,:)-E(i,:));
end
EST=mean(Y);
SE=std(Y);
HW=SE*icdf('normal',[0.975],0,1);
end