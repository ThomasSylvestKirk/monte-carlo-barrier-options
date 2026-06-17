function [EST,SE,HW,b] = CVforDOCCONT(S0,D,K,sigma,T,n,repetetions,r,Z)
%This functions simulates the option and use regress to estimate b to use control
%variates to estimate the payoff.
h=T/n;
repetetionsforest=repetetions/10;
DC=D*exp(-0.5826*sigma*sqrt(h));
[Y,X]=DOCCV(S0,D,K,sigma,T,n,repetetionsforest,r,Z(1:n,1:repetetionsforest),DC,h);
E=[1;DOCbls(S0,T,K,r,sigma,DC)];
b=regress(Y',X');
clear Y X
repetetions=repetetions-repetetionsforest;
[Y,X]=DOCCV(S0,D,K,sigma,T,n,repetetions,r,Z(1:n,1+repetetionsforest:end),DC,h);
for i=1+1:length(b)
Y=Y-b(i)*(X(i,:)-E(i,:));
end
EST=mean(Y);
SE=std(Y);
HW=SE*icdf('normal',[0.975],0,1);
end

function [Y,X]=DOCCV(S0,D,K,sigma,T,n,repetetions,r,Z,DC,h)
%ONLY USES CONTINUOUS VALUE. This function computes X and Y which have to b
%e used for the control variates
logS=nan(n+1,repetetions);
logS(1,:)=log(S0);
Y=nan(1,repetetions);
X=nan(1,repetetions);
LDC=log(DC);
LD=log(D);
for j=1:repetetions
    p=1; hit=0;
    for i=1:n
        logS(i+1,j)=logS(i,j)+h*(r-sigma^2/2)+sigma*Z(i,j)*sqrt(h);
            if logS(i+1,j)<LDC
                     Y(j)=0; X(j)=0; hit=1;
                     break
            elseif logS(i+1,j)<LD
                Y(j)=0;
                p=CSP(logS(:,j),i,LDC,sigma,h);
                     if i==n
                         X(j)=p*exp(-r*T)*(max(exp(logS(i+1,j))-K,0));
                     else
                         X(j)=p*exp(-r*i*h)*DOCbls(exp(logS(i+1,j)),(n-i)*h,K,r,sigma,DC);
                     end
                     hit=1; 
                     break
            end
    end
                if hit==0;
                    Y(j)=exp(-r*T)*max(exp(logS(end,j))-K,0);
                    p=CSP(logS(:,j),n,LDC,sigma,h);
                    X(j)=p*Y(j);
                end
end
X=[ones(1,repetetions);X];
end

function p=CSP(X,n,D,sigma,monitoringinterval)
%TRæk IID, D, and X should already be in log terms
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