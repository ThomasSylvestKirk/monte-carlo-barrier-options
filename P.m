function prob=P(X0,X,delta,D,sigma,monitoringinterval)
%This is used for the BB. Uses the brownian bridge tocalculate the trigger
%probability
prob=1;
    if X<D | X0<D
    else
        value=-2*((X0-D)*(X-D))/(sigma^2);
    if value>-36.044
        prob=exp(value);
    else 
    end
    end