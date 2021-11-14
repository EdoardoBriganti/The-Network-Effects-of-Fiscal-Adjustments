function [Bayesian_mcmc_results,Bayesian_mcmc_results_ind] = ave_eff_MC(rho,tau_or_gamma,...
    n,W,tax_or_exp,data,uncertainty,share,industry_weight_string,industry_weight,periods)


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Preliminary Objects:

% Calculate statistics about fiscal shocks (necessary for fiscal plans
% simulations)
[fisc_cons_ols,fisc_cons_summary] = fiscal_shock_analysis(data);

% Setting sizes:
k  = size(rho,1);
MC = 5*1e4;

% Determine how to weigh each industry when averaging the effects:
switch industry_weight_string
    case "average"
        weight = ones(1,n) ./ n;  % (1 x n)
    case "industry_share"
        weight = industry_weight; % (1 x n)
end
% Calculate the results by aggregating industry values:
[i,j] = size(weight);
if i>j
    weight = transpose(weight);
end



%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% MonteCarlo simulation of Average Effects:

% Preallocation
Bayesian_mcmc_results_ind = zeros(n,MC,3);
Bayesian_mcmc_results     = zeros(MC,3);

% MonteCarlo:
for mc = 1:MC
    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % Draw an iteration:
    draw = randi(k);
    % Draw fiscal coefficients from the posterior distributions
    mc_theta = transpose(tau_or_gamma(draw,:)); 
    % Draw spatial autoregressive parameter from the posterior distribution:
    mc_rho  = rho(draw);
    
    % Simulate Fiscal Plan:
    switch tax_or_exp
        case "tax"
            [mc_shock,~] = fiscal_shock_generator(fisc_cons_ols,fisc_cons_summary,uncertainty);
        case "exp"
            [~,mc_shock] = fiscal_shock_generator(fisc_cons_ols,fisc_cons_summary,uncertainty);
    end
    
    

    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=      
    % Calculate Results Industry by Industry:
    % Total:
    Total_ind1 = (eye(n)-mc_rho*W) \ share * (...
        mc_shock(1) * mc_theta(1) + sum(mc_shock(3:end)) * mc_theta(3)) ; 
    % Instantaneous:
    Instant_ind1 =  share * (...
        mc_shock(1) * mc_theta(1) + sum(mc_shock(3:end)) * mc_theta(3) );
    % Network:
    Network_ind1 = Total_ind1 - Instant_ind1;

    
    % Calculate Results by averaging across the industries:
    % Total
    Total1 = squeeze(weight * Total_ind1);
    % Direct:
    Instant1 = squeeze(weight * Instant_ind1);
    % Indirect:
    Network1 = squeeze(weight * Network_ind1);
    
    switch periods
        
        case "1year"
            
            %-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=
            % Calculate the effect over two years effect: overall percent
            % change:
            Total_ind   = Total_ind1;
            Instant_ind = Instant_ind1;
            Network_ind = Network_ind1;
            
            Total   = Total1;
            Instant = Instant1;
            Network = Network1;
            
            
        case "2years"
            
            %-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=
            % Calculate the effect of the implementation of the announced
            % fiscal plan, using the anticipated coefficient:
            
            % Total:
            Total_ind2 = (eye(n)-mc_rho*W) \ share * (sum(mc_shock(3:end)) * mc_theta(2)) ; 
            % Instantaneous:
            Instant_ind2 =  share * (sum(mc_shock(3:end)) * mc_theta(2) );
            % Network:
            Network_ind2 = Total_ind2 - Instant_ind2;

            % Calculate Results by averaging across the industries:
            % Total
            Total2 = squeeze(weight * Total_ind2);
            % Direct:
            Instant2 = squeeze(weight * Instant_ind2);
            % Indirect:
            Network2 = squeeze(weight * Network_ind2);
            
            %-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=
            % Calculate the effect over two years effect: overall percent
            % change:
            Total_ind   = Total_ind1 + Total_ind2;
            Instant_ind = Instant_ind1 + Instant_ind2;
            Network_ind = Network_ind1 + Network_ind2;
            
            Total   = Total1 + Total2;
            Instant = Instant1 + Instant2;
            Network = Network1 + Network2;
            
    end

    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % Store Results:
    Bayesian_mcmc_results_ind(:,mc,:) = cat(3,Total_ind,Instant_ind,Network_ind);
    Bayesian_mcmc_results(mc,:) = [Total Instant Network];

end

end