function [EST,SE,HW,b,Rsquared] = CVEURUND(S0,D,K,sigma,T,n,repetetions,r,Z)
%This functions simulates the option and use regress to estimate b. Uses control
%variates to estimate the payoff.
h=T/n;
repetetionsforest=repetetions/10;
[Y,X]=DOCEURUND(S0,D,K,sigma,T,n,repetetionsforest,r,Z(:,1:repetetionsforest));
E= [1;blsprice(S0,K,r,T,sigma);S0];
[b,~,~,~,stats]=regress(Y',X');
Rsquared=stats(1);
clear Y X
repetetions=repetetions-repetetionsforest;
[Y,X]=DOCEURUND(S0,D,K,sigma,T,n,repetetions,r,Z(:,repetetionsforest+1:end));

for i=1+1:length(b)
Y=Y-b(i)*(X(i,:)-E(i,:));
end
EST=mean(Y);
SE=std(Y);
HW=SE*icdf('normal',[0.975],0,1);
end