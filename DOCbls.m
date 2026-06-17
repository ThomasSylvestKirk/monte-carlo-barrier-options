function [DOC] = DOCbls(X0,T,K,r,sigma,D)
%     d1=(log(X0/K)+(r+(sigma^2)/2)*T)/(sigma*sqrt(T));
%     d2=d1-sigma*sqrt(T);
%     c=X0*normcdf(d1)-K*exp(-r*T)*normcdf(d2);
%     d11=(log((D^2/X0)/K)+(r+(sigma^2)/2)*T)/(sigma*sqrt(T));
%     d22=d11-sigma*sqrt(T);
%     C2=(D^2/X0)*normcdf(d11)-K*exp(-r*T)*normcdf(d22);    
%     DOC=c-(X0/D)^(1-2*r/(sigma^2))*C2;
%use the formulas from section 25.8 Hull (2012)    

if D<X0
lambda=(r)/(sigma^2)+1/2;
 x1=(log(X0/D))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
 y1=(log(D/X0))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
 DOC=X0*normcdf(x1)-K*exp(-r*T)*normcdf(x1-sigma*sqrt(T))-X0*(D/X0)^(2*lambda)*normcdf(y1)+K*exp(-r*T)*(D/X0)^(2*(lambda-1))*normcdf(y1-sigma*sqrt(T));
else
    lambda=(r)/(sigma^2)+1/2;
    y=log((D^2)/(X0*K))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
    d1=(log(X0/K)+(r+(sigma^2)/2)*T)/(sigma*sqrt(T));
    d2=d1-sigma*sqrt(T);
    c=X0*normcdf(d1)-K*exp(-r*T)*normcdf(d2);
    DOC=c-X0*(D/X0)^(2*lambda)*normcdf(y)-K*exp(-r*T)*(D/X0)^(2*(lambda-1))*normcdf(y-sigma*sqrt(T));
end


% if X0>D
% if D>=K
% lambda=(r)/(sigma^2)+1/2;
%  x1=(log(X0/D))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
%  y1=(log(D/X0))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
%  DOC=X0*normcdf(x1)-K*exp(-r*T)*normcdf(x1-sigma*sqrt(T))-X0*(D/X0)^(2*lambda)*normcdf(y1)+K*exp(-r*T)*(D/X0)^(2*(lambda-1))*normcdf(y1-sigma*sqrt(T));
% else
%     lambda=(r)/(sigma^2)+1/2;
%     y=log((D^2)/(X0*K))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
%     d1=(log(X0/K)+(r+(sigma^2)/2)*T)/(sigma*sqrt(T));
%     d2=d1-sigma*sqrt(T);
%     c=X0*normcdf(d1)-K*exp(-r*T)*normcdf(d2);
%     DOC=c-X0*(D/X0)^(2*lambda)*normcdf(y)-K*exp(-r*T)*(D/X0)^(2*(lambda-1))*normcdf(y-sigma*sqrt(T));
% end
% else
%     DOC=0;
end