function [rho_new,c_prop_new,acceptance_rate_new] = metropolis_hastings(rho,...
    c_prop_old,acceptance_rate_old,n,W,t_i,sig2,dy_i,Inv_V,beta,Xi,d,mc,a,b)


% vi: (n x 1) vector  of vi (drawn sector specific variances)

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Step 1) Update value of rho using Metroplis Hastings Algorithm:

% Construct Candidate Rho* using a random walk with variance equal to: c_prop^2.
rho_star = rho + c_prop_old*randn;

%  Discard the candidate if it does not belong to the support (a,b)
if rho_star<a || rho_star>b  
    
    % If outside of the support, then the new value is the old one:
    rho_new = rho;
    % Therefore we reject:
    s_i = 0;
    
else
    
    % Compute necessary objects for Psi:
    A_star = eye(n) - rho_star * W;
    A      = eye(n) - rho      * W;
    prior = ( ( (rho_star-a)*(b-rho_star) )/( (rho-a)*(b-rho) ) )^(d-1);
    SUM = 0;
    for t = 1:t_i
        yMy   = dy_i(:,t)' * ( A_star' * Inv_V * A_star - A' * Inv_V * A ) * dy_i(:,t);
        bxVMy = beta' *  Xi(:,:,t)' * Inv_V * (A_star - A) * dy_i(:,t);
        SUM = SUM + (yMy - 2*bxVMy) ;
    end
    
    % Compute Psi:
    Psi   = (det(A_star)/det(A))^(t_i) * exp( -1/(2*sig2) * SUM ) * prior;
    
    % Compute the probability of acceptance
    pi_i = min([1 Psi]);
    
    % Bernoulli experiment: accept or reject?
    s_i = binornd(1,pi_i);
    
    % Update the value of rho:
    rho_new = s_i * rho_star + (1-s_i) * rho;
end




%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Step 2) Update variance of the random walk: 

% Update Acceptance Rate:
acceptance_rate_new = ( acceptance_rate_old*(mc-1) + s_i )/mc;

% Update variance:
if acceptance_rate_new < 0.4 % Reduce variance if the acceptance rate is too low
    
    c_prop_new = c_prop_old/1.1; 
    
elseif acceptance_rate_new > 0.6 % Increase variance if the acceptance rate is too high
    
    c_prop_new = 1.1*c_prop_old; 
    
else
    c_prop_new = c_prop_old;
end
    
    
end