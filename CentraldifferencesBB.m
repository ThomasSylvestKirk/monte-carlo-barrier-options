function[sensitivity,varianceofsensitivity]=CentraldifferencesBB(S0,bumpS,bumpV,D,K,sigma,T,n,dim,r,Z)
%Finite difference method central difference
%This function computes the sensitivities of the barrier option using the
%Brownian bridge to estimate the options price to estimate the
%sensitivities
%Use a bump to one of the underlying parameters
sensitivity=0;
varianceofsensitivity=0;
bumps=0;
    [~,~,~,fm] = DOCBB(S0-bumpS,D,K,sigma-bumpV,T,n,dim,r,Z(1:n,:)); 
    [~,~,~,fp] = DOCBB(S0+bumpS,D,K,sigma+bumpV,T,n,dim,r,Z(1:n,:));
    %Estimate by central Finite diff
    Del = (fp-fm)/(2*(bumpS+bumpV)); %Get FD's
    val =  sum(Del)/dim;      %Estimate by average
    var = (sum(Del.^2)/dim - val^2)/(dim-1); %Variance estimate
    %Keeps track of estimates, variances and bump sizes
    sensitivity  = [val];
    varianceofsensitivity  = [var];
end
