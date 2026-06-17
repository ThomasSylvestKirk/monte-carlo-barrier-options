function [Y,X] = DOC3CV(S0,D,K,sigma,T,n,repetetions,r,Z)
%Uses thre control variates. This function computes X and Y which have to b
%e used for the control variates
logS=nan(n+1,repetetions);
logS(1,:)=log(S0);
h=T/n;
LD=log(D);
Y=nan(1,repetetions);
X=nan(2,repetetions);

for j=1:repetetions
    p=1; hit=0;
    for i=1:n
        logS(i+1,j)=logS(i,j)+h*(r-sigma^2/2)+sigma*Z(i,j)*sqrt(h);
            if logS(i+1,j)<LD
                     Y(j)=0;  hit=1;            
            end
    end
                if hit==0;
                    Y(j)=exp(-r*T)*max(exp(logS(end,j))-K,0);
                end
end
S(1:n,1:repetetions) = cumsum((r-0.5*(sigma.^2))*h+sigma*sqrt(h).*Z);
S(:,:) = S+logS(1,:);
S=exp(S);
X(1,:)=exp(-r*T)*max(S(end,:)-K,0);
X(2,:)=exp(-r*T)*S(end,:); 
X=[ones(1,repetetions);X];
end