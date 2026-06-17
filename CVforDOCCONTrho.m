function [EST,SE,HW,Y,b,rho] = CVforDOCCONTrho(S0,D,K,sigma,T,n,repetetions,r,Z)
%This functions simulates the option and use regress to estimate b and rho
%of the control variates. It is splitted into a separate function to
%minimize the TF 
h=T/n;
repetetionsforest=repetetions/10;
[Y,X]=DOCCVrho(S0,D,K,sigma,T,n,repetetionsforest,r,Z(1:n,1:repetetionsforest));
DC=D*exp(-0.5826*sigma*sqrt(h));
E=[1;DOCbls(S0,T,K,r,sigma,DC)];
b=regress(Y',X');
clear Y X
repetetions=repetetions-repetetionsforest;
[Y,X]=DOCCVrho(S0,D,K,sigma,T,n,repetetions,r,Z(1:n,1+repetetionsforest:end));
VARY=var(Y);
for i=1+1:length(b)
Y=Y-b(i)*(X(i,:)-E(i,:));
end
VARYY=var(Y);
rho=1-VARYY/VARY;
EST=mean(Y);
SE=std(Y);
HW=SE*icdf('normal',[0.975],0,1);
end