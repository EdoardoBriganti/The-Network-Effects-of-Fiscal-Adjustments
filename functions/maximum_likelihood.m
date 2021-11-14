function [MLE,rho_down,rho_up,beta,Omega] = maximum_likelihood(data_MLE,model,inversion)


% Extract Data for MLE Estimation:
T       = data_MLE.T;
n       = data_MLE.n;
dy      = data_MLE.dy;
d_down  = data_MLE.d_down;
d_up    = data_MLE.d_up;
Adown   = data_MLE.Adown;
Aup     = data_MLE.Aup;
X       = data_MLE.X;
lb_up   = data_MLE.lb_up;
ub_up   = data_MLE.ub_up;
lb_down = data_MLE.lb_down;
ub_down = data_MLE.ub_down;


% Determine number of regressors per equation except for the spatial variable:
switch model
    case "static"
        k_eq = 1 + 6 + 1;
    case "dynamic"
        k_eq = 1 + 1 + 6 + 1;
end


% Construct Log_Likelihood:
Log_L   = @(r) log_likelihood(r,T,n,dy,d_down,d_up,Adown,Aup,X,k_eq);  
% Don't show messages while optimizing:
options = optimoptions('fmincon','display','off');
% Initial value
x0      = .1 + .9*rand(1,2); 
% Boundaries of spatial coeffiecients:
lb      = [lb_down,lb_up];
ub      = [ub_down,ub_up]; 
% Constrained Optimization:
sol = fmincon(Log_L,x0,[],[],[],[],lb,ub,[],options);

% Calculate the other parameters: beta and Omega
[~,beta,Omega] = log_likelihood(sol,T,n,dy,d_down,d_up,Adown,Aup,X,k_eq);

% Determine Indices of hat_beta to reconnect them to their parameters. This
% order is established by the order of appearance of the variables in the
% regressor matrix X:
index_FE    = transpose(1:n);
index_tau   = transpose(index_FE(end)+1:index_FE(end)+3);
index_gamma = transpose(index_tau(end)+1:index_tau(end)+3);
index_08_09 = transpose(index_gamma(end)+1:index_gamma(end)+2);
switch model
    case "static"
        index_phi   = [];
    case "dynamic"
        index_phi   = transpose(index_08_09(end)+1:index_08_09(end)+n);
end

% Rename Solutions:
rho_down = sol(1);
rho_up   = sol(2);

% Storing results:
switch inversion
    case "no"
        MLE  = [rho_down ; beta(index_tau) ; ...
            rho_up ; beta(index_gamma);...
            beta(index_08_09);...
            beta(index_FE) ; beta(index_phi); diag(Omega)];
    case "yes"
        MLE  = [rho_up ; beta(index_tau) ; ...
            rho_down ; beta(index_gamma);...
            beta(index_08_09);...
            beta(index_FE) ; beta(index_phi); diag(Omega)];
end


% Calculate Pointwise Log_Likelihood: l_t
eps_t = (eye(n) - rho_down*Adown*d_down - rho_up*Aup*d_up) * dy_m - X_m * beta;
Log_L_t  = -T/2 * log(det(omega)) + t_down * log(det(eye(n)-rho_down*Adown)) + ...
            + t_up * log(det(eye(n)-rho_up*Aup)) ;

end