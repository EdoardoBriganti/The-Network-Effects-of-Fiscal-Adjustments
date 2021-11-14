function tab = table_par_estimation(posterior,share_EB,share_TB,n,varnames,inversion,model,tab_data_MLE)



% Number of variables:
nVar = size(varnames,1);

% Index of variables in beta:
index_FE    = (1:n);
index_tau   = (index_FE(end)+1:index_FE(end)+3);
index_gamma = (index_tau(end)+1:index_tau(end)+3);
index_08_09 = (index_gamma(end)+1:index_gamma(end)+2);
switch model
    case "static"
        index_phi = [];
    case "dynamic"
        index_phi = (index_08_09(end)+1:index_08_09(end)+n);
end

% Store results of Bayesian MCMC:
switch inversion
    case "no"
        tab_data = [ posterior.rho_down ...
            posterior.beta(:,index_tau) .* mean(share_TB) ...
            posterior.rho_up ...
            posterior.beta(:,index_gamma) .* mean(share_EB) ...
            posterior.beta(:,index_08_09) ...
            posterior.beta(:,index_phi) posterior.beta(:,index_FE) ...
            posterior.sig2 .* transpose(posterior.vi)];
    case "yes"
        tab_data = [ posterior.rho_up ...
            posterior.beta(:,index_tau) .* mean(share_TB) ...
            posterior.rho_down ...
            posterior.beta(:,index_gamma) .* mean(share_EB) ...
            posterior.beta(:,index_08_09) ...
            posterior.beta(:,index_phi) posterior.beta(:,index_FE) ...
            posterior.sig2 .* transpose(posterior.vi)];
end


% Iteration of Bayesian MCMC:
nIter = size(tab_data,1);


% Adjust fiscal coefficients to the mean-share:
% TB:
tab_data_MLE(2:4,:) = tab_data_MLE(2:4,:) .* mean(share_TB);
% EB:
tab_data_MLE(6:8,:) = tab_data_MLE(6:8,:) .* mean(share_EB);


% Construct the table:
tab_m = zeros(nVar,12);
for j = 1 : nVar
    
    % MLE Point Estimate and Std constructed using the Analytical Fisher Information:
    tab_m(j,1:2) = tab_data_MLE(j,:);
    % Mean of the Posterior:
    tab_m(j,3) = mean(tab_data(:,j));
    % Std of the Posterior:
    tab_m(j,4) = std(tab_data(:,j));
    % Pr(theta<0|D):
    tab_m(j,5) = sum(tab_data(:,j)<0)/nIter;
    % Quantiles of the posterior:
    tab_data(:,j) = sort(tab_data(:,j));
    tab_m(j,6) = tab_data(round(0.05*nIter),j);
    tab_m(j,7) = tab_data(round(0.1*nIter),j);
    tab_m(j,8) = tab_data(round(0.16*nIter),j);
    tab_m(j,9) = tab_data(round(0.5*nIter),j);
    tab_m(j,10) = tab_data(round(0.84*nIter),j);
    tab_m(j,11) = tab_data(round(0.9*nIter),j);
    tab_m(j,12) = tab_data(round(0.95*nIter),j);
end
tab = array2table(tab_m,'VariableNames',{'MLE','MLE std','Mean','std',...
    'prob_neg','q5','q10','q16','q50','q84','q90','q95',});
tab = addvars(tab,varnames,'Before','MLE');


end