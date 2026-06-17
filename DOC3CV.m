function [Y,X] = DOC3CV(S0,D,K,sigma,T,n,repetetions,r,Z)
%Uses thre control variates. This function computes X and Y which have to b
%e used for the control variates
logS=nan(n+1,repetetions);
logS(1,:)=log(S0);
h=T/n;
DC=D*exp(-0.5826*sigma*sqrt(h));
LDC=log(DC);
LD=log(D);
Y=nan(1,repetetions);
X=nan(3,repetetions);

for j=1:repetetions
    p=1; hit=0;
    for i=1:n
        logS(i+1,j)=logS(i,j)+h*(r-sigma^2/2)+sigma*Z(i,j)*sqrt(h);
            if logS(i+1,j)<LDC
                     Y(j)=0; X(1,j)=0; hit=1;
            elseif logS(i+1,j)<LD
                Y(j)=0;
                p=CSP(logS(:,j),i,LDC,sigma,h);
                     if i==n
                         X(1,j)=p*exp(-r*T)*(max(exp(logS(end,j))-K,0)); 
                     else
                         X(1,j)=p*exp(-r*i*h)*DOCbls(exp(logS(i+1,j)),(n-i)*h,K,r,sigma,DC);                       
                     end
                     hit=1;               
            end
    end
                if hit==0;
                    Y(j)=exp(-r*T)*max(exp(logS(end,j))-K,0);
                    p=CSP(logS(:,j),n,LDC,sigma,h);
                    X(1,j)=p*Y(j);
                end
end
% S(1:n,1:repetetions) = cumsum((r-0.5*(sigma.^2))*h+sigma*sqrt(h).*Z);
% S(:,:) = S+log(S0);
S=exp(logS(end,:));
X(2,:)=exp(-r*T)*max(S(end,:)-K,0);
X(3,:)=exp(-r*T)*S(end,:);  
X=[ones(1,repetetions);X];
end