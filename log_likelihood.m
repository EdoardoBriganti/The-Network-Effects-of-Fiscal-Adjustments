function [Log_L,beta_wls,omega] = log_likelihood(r,T,n,dy,d_down,d_up,Adown,Aup,X,k_eq)

%{
        PURPOSE: Calculate the log likelihood of the model.

%}

% NUmber of years of fiscal consolidations:
t_down = sum(d_down);
t_up = sum(d_up);

% Reshape Dummies:
D_down = kron(d_down,ones(n,1));
D_up = kron(d_up,ones(n,1));


% Concentrate the Log_Likelihood:
Z        = dy - r(1)*D_down.*(kron(eye(T),Adown)*dy) - r(2)*D_up.*(kron(eye(T),Aup)*dy);
beta_ols = X\Z;
eps_ols  = Z - X*beta_ols; 
omega    = diag( 1/(T-k_eq).*sum(reshape(eps_ols,[n,T]).*reshape(eps_ols,[n,T]),2) );
sigma    = kron(eye(T),omega);
beta_wls = ( X'*(sigma\X) )\( X'*( sigma\Z) );
eps_wls  = Z - X*beta_wls;
Hdown       = inv(eye(n) - r(1)*Adown);
Hup       = inv(eye(n) - r(2)*Aup);

% Calculate Log_likelihood:
Log_L    = -( - T/2*log(det(omega)) - t_down*log(det(Hdown)) ...
    - t_up*log(det(Hup)) - 0.5*eps_wls'*(sigma\eps_wls) ) ;

end