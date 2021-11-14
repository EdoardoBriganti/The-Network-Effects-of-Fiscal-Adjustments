function posterior = cum_irf_MC(posterior,n,MC,horz,Adown,Aup,share_TB,share_EB,...
    data,uncertainty,inversion)

% Derive descriptive statistics of original fiscal consolidations, which
% will be employed to simulate "reasonable" fiscal plans, in the sense that
% they will mimic the average in-sample fiscal plan:
[b,summary,~] = fiscal_shock_analysis(data);


% Index of variables in beta:
index_FE    = (1:n);
index_tau   = (index_FE(end)+1:index_FE(end)+3);
index_gamma = (index_tau(end)+1:index_tau(end)+3);
index_08_09 = (index_gamma(end)+1:index_gamma(end)+2);
index_phi   = transpose(index_08_09(end)+1:index_08_09(end)+n);

% Extract the data:
rho_down = posterior.rho_down; 
rho_up   = posterior.rho_up;
tau      = posterior.beta(:,index_tau);
gamma    = posterior.beta(:,index_gamma);
phi      = posterior.beta(:,index_phi);


% Number of elements in the posteriors:
nMCMC = size(rho_down,1);


% Preallocate Estimation Target:
% TB:
cum_irf_tax_tot = zeros(n,horz+1,MC);
cum_irf_tax_dir = zeros(n,horz+1,MC);
cum_irf_tax_net = zeros(n,horz+1,MC);
% EB:
cum_irf_exp_tot = zeros(n,horz+1,MC);
cum_irf_exp_dir = zeros(n,horz+1,MC);
cum_irf_exp_net = zeros(n,horz+1,MC);


% Running an MC simulation to construct the posterior distributions of the 
% "dynamic Average Effects":
parfor mc = 1 : MC
    
    
    % Step 1) draw parameters from the posteriors. Draw the parameters from
    %         the same iteration in the Bayesian MCMC. 
    draw = randi(nMCMC);
    
    mc_rho_down = rho_down(draw);
    mc_rho_up   = rho_up(draw);
    mc_Phi    = diag(phi(draw,:));
    mc_tau    = tau(draw,:);
    mc_gamma  = gamma(draw,:);
    
    
    % Step 2) Simulate Fiscal Plan:
    [tax,exp] = fiscal_shock_generator(b,summary,uncertainty);
    % Collapse the future three years shocks into a single year future
    % shock:
    tax = [tax(1:2) ; sum(tax(3:5))];
    exp = [exp(1:2) ; sum(exp(3:5))];
    
    % Step 3) Constructing IRF:
    shock_tax = zeros(horz+1,1);
    shock_exp = zeros(horz+1,1);
    
    irf_tax_tot = zeros(n,horz+1);
    irf_exp_tot = zeros(n,horz+1);
    irf_tax_dir = zeros(n,horz+1);
    irf_exp_dir = zeros(n,horz+1);
    irf_tax_net = zeros(n,horz+1);
    irf_exp_net = zeros(n,horz+1);

    for t = 1:horz
        if t == 1
            % Tax
            shock_tax(t+1) = mc_tau(1)*tax(1) + mc_tau(2)*tax(2) + mc_tau(3)*tax(3);
            % Exp
            shock_exp(t+1) = mc_gamma(1)*exp(1) + mc_gamma(2)*exp(2) + mc_gamma(3)*exp(3);
            
        elseif t == 2
            % Tax 
            shock_tax(t+1) = mc_tau(2)*tax(3);
            % Exp
            shock_exp(t+1) = mc_gamma(2)*exp(3);
  
        end
        
        
        % Calculate IRF: 
        switch inversion
            case "no"
                M_tau   = (eye(n) - mc_rho_down * Adown);
                M_gamma = (eye(n) - mc_rho_up   * Aup);
                
            case "yes"
                M_tau   = (eye(n) - mc_rho_up   * Aup);
                M_gamma = (eye(n) - mc_rho_down * Adown);
                
        end
        
        % Tax IRF:
        irf_tax_tot(:,t+1) = M_tau \ ...
            ( mc_Phi * irf_tax_tot(:,t) + share_TB * shock_tax(t+1) );
        irf_tax_dir(:,t+1) = mc_Phi * irf_tax_dir(:,t) + share_TB * shock_tax(t+1);
        irf_tax_net(:,t+1) = irf_tax_tot(:,t+1) - irf_tax_dir(:,t+1);
        
        irf_exp_tot(:,t+1) = M_gamma \ ...
            ( mc_Phi * irf_exp_tot(:,t) + share_EB * shock_exp(t+1) );
        irf_exp_dir(:,t+1) = mc_Phi * irf_exp_dir(:,t) + share_EB * shock_exp(t+1);
        irf_exp_net(:,t+1) = irf_exp_tot(:,t+1) - irf_exp_dir(:,t+1);
        
    end
    
    % Step 5: Calculate Cumulative IRFs.
    % TB Plans:
    cum_irf_tax_tot(:,:,mc) = cumsum(irf_tax_tot,2);
    cum_irf_tax_dir(:,:,mc) = cumsum(irf_tax_dir,2);
    cum_irf_tax_net(:,:,mc) = cumsum(irf_tax_net,2);
    
    % TB Plans:
    cum_irf_exp_tot(:,:,mc) = cumsum(irf_exp_tot,2);
    cum_irf_exp_dir(:,:,mc) = cumsum(irf_exp_dir,2);
    cum_irf_exp_net(:,:,mc) = cumsum(irf_exp_net,2);

    
end


% Store results
posterior.cum_irf_tax_tot = cum_irf_tax_tot;
posterior.cum_irf_tax_dir = cum_irf_tax_dir;
posterior.cum_irf_tax_net = cum_irf_tax_net;
posterior.cum_irf_exp_tot = cum_irf_exp_tot;
posterior.cum_irf_exp_dir = cum_irf_exp_dir;
posterior.cum_irf_exp_net = cum_irf_exp_net;

end