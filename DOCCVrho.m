function [Y,X]=DOCCV(S0,D,K,sigma,T,n,repetetions,r,Z)
%ONLY USES CONTINUOUS VALUE. This function computes X and Y which are used to calculate b
logS=nan(n+1,repetetions);
logS(1,:)=log(S0);
h=T/n;
DC=D*exp(-0.5826*sigma*sqrt(h));
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
                         X(j)=p*exp(-r*(i)*h)*DOCbls(exp(logS(i+1,j)),(n-i)*h,K,r,sigma,DC); % Conditional Barrierprice opret fil
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