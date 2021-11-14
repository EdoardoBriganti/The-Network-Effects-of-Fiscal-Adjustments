function [MLE_std,MLE_varcov,FIm] = fisher_information_matrix(data_MLE,k,...
    X_down,X_up,rho_down,rho_up,beta,Omega,inversion,model)
%{
        PURPOSE: 
        Calculate the Fisher Information Matrix. 

%}

% Extract Data:
T       = data_MLE.T;
n       = data_MLE.n;
d_down  = data_MLE.d_down;
d_up    = data_MLE.d_up;
Adown   = data_MLE.Adown;
Aup     = data_MLE.Aup;
X       = data_MLE.X;

% Calculate Years of fiscal consolifations, that is, years when the
% downstream and upstream network are activated:
t_down = sum(d_down);
t_up   = sum(d_up);

% Constructing useful elements:
A1      = inv(eye(n)-rho_down*Adown);
A2      = inv(eye(n)-rho_up*Aup);

M1      = A1'*Adown'*inv(Omega)*Adown*A1;
M2      = inv(Omega)*Adown*A1;
M3      = A2'*Aup'*inv(Omega)*Aup*A2;
M4      = inv(Omega)*Aup*A2;




%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%  Calculating the FIM: 

FIm            = zeros(2+k+n,2+k+n);
% First row:
FIm(1,1)       = t_down*trace(Adown*A1*Adown*A1 + Omega*M1) + beta'*X_down'*kron(eye(t_down),M1)*X_down*beta;
FIm(1,2)       = 0;
FIm(1,3:k+2)   = beta'*X_down'*kron(eye(t_down),M2')*X_down;
FIm(1,k+3:end) = t_down*diag(M2)';
% Second row:
FIm(2,1)       = 0;
FIm(2,2)       = t_up*trace(Aup*A2*Aup*A2 + Omega*M3) + beta'*X_up'*kron(eye(t_up),M3)*X_up*beta;
FIm(2,3:k+2)   = beta'*X_up'*kron(eye(t_up),M4')*X_up;
FIm(2,k+3:end) = t_up*diag(M4)';
% Third to 23rd row
FIm(3:k+2,1)       = FIm(1,3:k+2)';
FIm(3:k+2,2)       = FIm(2,3:k+2)';
FIm(3:k+2,3:k+2)   = X'*kron(eye(T),inv(Omega))*X;
FIm(3:k+2,k+3:end) = zeros(k,n);
% 24th to 38th row (end)
FIm(k+3:end,1)       = FIm(1,k+3:end)';
FIm(k+3:end,2)       = FIm(2,k+3:end)';
FIm(k+3:end,3:k+2)   = FIm(3:k+2,k+3:end)';
FIm(k+3:end,k+3:end) = T/2 * diag(diag(inv(Omega)).^2);



%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Computing the covariance matrix of the MLE (inverting FIM):
MLE_varcov  = inv(FIm);


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Computing standard deviations of the estimated parameters:
MLE_std = sqrt(diag(MLE_varcov));


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Reorder MLE_std to make it comparable with the order of the MLE:
% Storing results:

% Determine Indices of hat_beta to reconnect them to their parameters. This
% order is established by the order of appearance of the variables in the
% regressor matrix X:
index_FE    = transpose(2+1:2+n);
index_tau   = transpose(index_FE(end)+1:index_FE(end)+3);
index_gamma = transpose(index_tau(end)+1:index_tau(end)+3);
index_08_09 = transpose(index_gamma(end)+1:index_gamma(end)+2);
switch model
    case "static"
        index_phi   = [];
        index_sig2  = transpose(index_08_09(end)+1:index_08_09(end)+n);
    case "dynamic"
        index_phi   = transpose(index_08_09(end)+1:index_08_09(end)+n);
        index_sig2  = transpose(index_phi(end)+1:index_phi(end)+n);
end

% Consruct the vector of MLE_Std according to the order of MLE:
switch inversion
    case "no"
        MLE_std  = [MLE_std(1) ; MLE_std(index_tau) ; ...
            MLE_std(2) ; MLE_std(index_gamma) ; ...
            MLE_std(index_08_09);...
            MLE_std(index_FE) ; MLE_std(index_phi) ; MLE_std(index_sig2)];
    case "yes"
        MLE_std  = [MLE_std(2) ; MLE_std(index_tau) ; ...
            MLE_std(1) ; MLE_std(index_gamma) ; ...
            MLE_std(index_08_09);...
            MLE_std(index_FE) ; MLE_std(index_phi) ; MLE_std(index_sig2)];
end

end