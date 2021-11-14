function tab = table_MLE(MLE,MLE_std,inversion,model,n,share_TB,share_EB)

% Construct string of variables' names:
FE   = [];
phi  = [];
sig2 = [];
for i = 1:n
    fe_str   = "alpha_" + string(i);
    phi_str  = "phi_"   + string(i);
    sig2_str = "sig2_"  + string(i);
    FE = [FE ; fe_str];
    switch model
        case "dynamic"
            phi = [phi ; phi_str];
    end
    sig2 = [sig2 ; sig2_str];
end

switch inversion
    case "no"
        variables = ["rho_down" ; "tau_u" ; "tau_a" ; "tau_f" ; ...
            "rho_up" ; "gamma_u" ; "gamma_a" ; "gamma_f" ; ...
            "D2008" ; "D2009" ; ...
            FE ; phi ; sig2];
    case "yes"
        variables = ["rho_up" ; "tau_u" ; "tau_a" ; "tau_f" ; ...
            "rho_down" ; "gamma_u" ; "gamma_a" ; "gamma_f" ; ...
            "D2008" ; "D2009" ; ...
            FE ; phi ; sig2];
end

% Adjust fiscal coefficients to the mean-weight:
% TB:
MLE(2:4)     = MLE(2:4) .* mean(share_TB);
MLE_std(2:4) = MLE_std(2:4) .* mean(share_TB);
% EB:
MLE(6:8)     = MLE(6:8) .* mean(share_EB);
MLE_std(6:8) = MLE_std(6:8) .* mean(share_EB);

% Construct Table:
tab = table(variables,MLE,MLE_std);

end