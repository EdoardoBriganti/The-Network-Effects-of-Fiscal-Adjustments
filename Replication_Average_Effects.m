%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   REPLICATION OF RESULTS - AVERAGE EFFECTS    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
PURPOSE: Replicate the results of the paper.




%=========================================================================%

Favero, Karamysheva, Briganti 
"The network effects of fiscal adjustments" (2021)

%=========================================================================%

%}


% Preamble:

clc
clear


% 1) Network matrices:
normalization = "no";    % "yes"    "no"
diagonal      = "with";  % "with"  "without"


% 2) Decide whether to invert the model, i.e. switch the dummies:
inversion     = "no";    % "yes"    "no"


% 3) Decide whether to run a dynamic or static simulation:
model = "static";   %  "static"     "dynamic"

% 4) Decides if you run the static model to simulate also the
% implementation of the future shock in year 2:
switch model
    case "static"
        % Decide how many periods to simulate: either 1 year or 2 years:
        periods = "2years";  % "1year"  "2years'
end

% 5) Determine whether to provide average effects or weighted average effects
% where the weights are represented by average industry shares conditional
% on TB and EB occurrence.
industry_weight_string = "industry_share";  % "average" "industry_share"
        
        
% 6) Decides wheter to run a placebo simulation or not:
placebo = "no";   %  "yes"    "no"
% Length of the Placebo:
nPl = 500;


% 7) Other:
% Fiscal Database:
cutoff = 50;
% Number of industries to include in the analysis (level of disaggregation)
num_ind = 62; %    "15"    "19"    "62"
% Fiscal database (shock distribution)
uncertainty   = 'yes';



%% Add path to Matlab Functions Subdirectory:
%addpath('/ULTIMATE CODE/Functions');



%% Upload variables of interest:


% Upload database:
[data,T,n,tb,eb,dy_m,A,A_hat_tr,share_EB,share_TB] = data_upload(cutoff,num_ind);


% IMPORTING GDP DATA: 
% dy_m is a (n x T) matrix - first year 1977, last year 2014
% Reorganize observatins into a vector (NT x 1)
dy       = reshape(dy_m,[T*n,1]);

% Constructing Lagged Value Added: matrix (n x T) - first year 1977, last year 2013
dy_lag_m = data.dy_real(:,2:end-1); 
% Constructing lags observations
dy_lag   = reshape(dy_lag_m,[T*n,1]);


% Upload and prepare the I-O matrices
[Adown,lb_down,ub_down]  = network_preparation(A,normalization,diagonal);
[Aup,lb_up,ub_up]  = network_preparation(A_hat_tr,normalization,diagonal);


% Controls for recessions to improve efficiency:
d_2008 = zeros(T,1);
d_2008(data.years==2008) = 1;
D_2008 = kron(d_2008,ones(n,1));

d_2009 = zeros(T,1);
d_2009(data.years==2009) = 1;
D_2009 = kron(d_2009,ones(n,1));

% Regressor matrix:
switch model
    case "static"
        X = [repmat(eye(n),T,1) ...
            kron(data.SU   .* tb ,share_TB) ...
            kron(data.SA   .* tb ,share_TB) ...
            kron(data.SA_F .* tb ,share_TB) ...
            kron(data.SU   .* eb ,share_EB) ...
            kron(data.SA   .* eb ,share_EB) ...
            kron(data.SA_F .* eb ,share_EB) ...
            D_2008 D_2009];
        
    case "dynamic"
        X = [repmat(eye(n),T,1) ...
            kron(data.SU   .* tb ,share_TB) ...
            kron(data.SA   .* tb ,share_TB) ...
            kron(data.SA_F .* tb ,share_TB) ...
            kron(data.SU   .* eb ,share_EB) ...
            kron(data.SA   .* eb ,share_EB) ...
            kron(data.SA_F .* eb ,share_EB) ...
            D_2008 D_2009 ...
            repmat(eye(n),[T,1]).*dy_lag ];
end


% Decide whether to run an inverted model: construct dummies D_down and
% D_up which are the dummies TB and EB which will be interacted with the
% downstream or upstream network.
switch inversion
    case "no"
        d_down = tb;
        d_up = eb;
    case "yes"
        d_down = eb;
        d_up = tb;
end
D_down = kron(d_down,ones(n,1));
D_up = kron(d_up,ones(n,1));

% Construct number of years in which we have TB and EB fiscal
% consolidations, which will be correspond to years when we have downstream
% or upstream spillovers:
t_down = sum(d_down);
t_up   = sum(d_up);


% Calculating regressors matrices in fiscal shocks years only:
k      = size(X,2);
XX     = permute(reshape(X',[k,n,T]),[2 1 3]);
X_down = transpose(reshape(permute(XX(:,:,find(d_down)),[2,1,3]),[k,n*t_down]));
X_up   = transpose(reshape(permute(XX(:,:,find(d_up)),[2,1,3]),[k,n*t_up]));


% Calculate same objects but in 3D matrix format with the third dimension
% represented by time:
% Reshape Regressors:
for t = 1:T
    
    % Initialize M:
    M = [];
    
    % Construct TB Fiscal Adjustments:
    TB_plans = [data.SU(t) data.SA(t) data.SA_F(t)] .* tb(t) .*share_TB;
    
    % Construct EB Fiscal Adjustments:
    EB_plans = [data.SU(t) data.SA(t) data.SA_F(t)] .* eb(t) .*share_EB;
    
    % Constgruct year dummies for the financial crisis:
    financial_crisis = [d_2008(t) d_2009(t)] .* ones(n,1);
    
    % Construct Regressor
    M = [eye(n) TB_plans EB_plans financial_crisis];
    
    % Include lag if model is dynamic:
    switch model
        case "dynamic"
            M = [M diag(dy_lag(:,t))]; 
    end
    
    % Store the 3D version of our dataset:
    X_m(:,:,t)      = M;
    X_down_m(:,:,t) = M .* d_down(t);
    X_up_m(:,:,t)   = M .* d_up(t);
end

% Remove "floors" which do not correspond to years of TB fiscal plans (EB
% if the model is inverted):
X_down_m(:,:,~any(X_down_m,[1,2])) = [];

% Remove "floors" which do not correspond to years of EB fiscal plans (TB
% if the model is inverted):
X_up_m(:,:,~any(X_up_m,[1,2])) = [];


%% Maximum Likelihood Estimation and Analytical Fisher Information Matrix:

% Prepare Data for MLE:
data_MLE.T = T;
data_MLE.n = n;
data_MLE.dy = dy;
data_MLE.d_down = d_down;
data_MLE.d_up = d_up;
data_MLE.Adown = Adown;
data_MLE.Aup = Aup;
data_MLE.X = X;
data_MLE.lb_up = lb_up;
data_MLE.ub_up = ub_up;
data_MLE.lb_down = lb_down;
data_MLE.ub_down = ub_down;

% Estimate parameters via ML:
[MLE,rho_down,rho_up,beta,Omega] = maximum_likelihood(data_MLE,model,inversion);
% Order of MLE depends on the inversion of the model

% Construct standard deviations using the analytical Fisher Information
% matrix:
MLE_std = fisher_information_matrix(data_MLE,k,X_down,X_up,rho_down,rho_up,beta,Omega,inversion,model);

% Summarize results in a Table:
table_ML = table_MLE(MLE,MLE_std,inversion,model,n,share_TB,share_EB);





%% Bayesian MCMC Estimation:

% Prepare data for Bayesian MCMC by saving all the necassary information
% into a cell array (cell arrays are better than structures since they
% allow for parallel loops:
data_MCMC = {dy_m ...
    d_down d_up t_down t_up ...
    k ...
    Adown Aup ub_up ub_down  ...
    X_m X_down_m X_up_m};


% Set length of final vector of Posterior Distribution:
nMCMC = 6*1e4;


switch placebo
    case "no"
        %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        % Estimate the parameters of the model via Bayesian MCMC:
        tic
        posterior = Bayesian_MCMC(data_MCMC,nMCMC,n,T);
        toc
        
        %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        % Construct a table which reports the estimation results, and 
        % compare them to the MLE ones.

        % Store name of variables from MLE table:
        varnames = table2array(table_ML(:,1));
        % Store results of MLE:
        tab_data_MLE = [MLE MLE_std];

        % Construct the Table:
        table_IIIa = table_par_estimation(posterior,share_EB,share_TB,n,...
            varnames,inversion,model,tab_data_MLE);

        
        
        
        %% Average Effects:
        
        % Indices within beta of fiscal coefficients:
        index_tau   = (n+1:n+3);
        index_gamma = (n+4:n+6);




        % Calculate Average Effects:
        switch model
            case "static"

                switch inversion
                    case "no"
                        % Simulation of average effects of a TB plan:
                        [posterior.AE_tax,posterior.AE_tax_ind] = ave_eff_MC(posterior.rho_down,...
                            posterior.beta(:,index_tau),n,Adown,"tax",data,uncertainty,...
                            share_TB,industry_weight_string,data.ind_weight_TB,periods);

                        % Simulation of average effects of an EB plan:
                        [posterior.AE_exp,posterior.AE_exp_ind] = ave_eff_MC(posterior.rho_up,...
                            posterior.beta(:,index_gamma),n,Aup,"exp",data,uncertainty,...
                            share_EB,industry_weight_string,data.ind_weight_EB,periods);

                    case "yes"

                        % Simulation of average effects of a TB plan:
                        [posterior.AE_tax,posterior.AE_tax_ind] = ave_eff_MC(posterior.rho_up,...
                            posterior.beta(:,index_tau),n,Aup,"tax",data,uncertainty,...
                            share_TB,industry_weight_string,data.ind_weight_TB,periods);

                        % Simulation of average effects of an EB plan:
                        [posterior.AE_exp,posterior.AE_exp_ind] = ave_eff_MC(posterior.rho_down,...
                            posterior.beta(:,index_gamma),n,Adown,"exp",data,uncertainty,...
                            share_EB,industry_weight_string,data.ind_weight_EB,periods);
                end

                %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                % Report results in a Table:
                tableIIIb = table_average_effect(posterior);

            case "dynamic"

                % Determine length of MonteCarlo:
                MC = 1e5;
                % Determine the number of years to simulate: horizon of the cumulative IRF:
                horz = 8;
                % Calculate cumulative impulse response:
                cum_irf = cum_irf_MC(posterior,n,MC,horz,Adown,Aup,share_TB,share_EB,data,uncertainty,inversion);


                %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                % Plot Cumulative Impulse Response Functions:
                fig = cumIRF(cum_irf);

        end
        
    case "yes"
        
        %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        % Step 1) Calculate results for original model:
        
        % Indices within beta of fiscal coefficients:
        index_tau   = (n+1:n+3);
        index_gamma = (n+4:n+6);

        
        %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        % Estimate parameters via Bayesian MCMC using original network:
        posterior = Bayesian_MCMC(data_MCMC,nMCMC,n,T);

        %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        % Calculate Average Effects using original network:
        
        % Simulation of average effects of a TB plan:
        AE_tax_orig = ave_eff_MC(posterior.rho_down,...
            posterior.beta(:,index_tau),n,Adown,"tax",data,uncertainty,...
            share_TB,industry_weight_string,data.ind_weight_TB,periods);
        % Store network Effect of TB Plans:
        AE_tax_net_orig = AE_tax_orig(:,3);

        % Simulation of average effects of an EB plan:
        AE_exp_orig = ave_eff_MC(posterior.rho_up,...
            posterior.beta(:,index_gamma),n,Aup,"exp",data,uncertainty,...
            share_EB,industry_weight_string,data.ind_weight_EB,periods);
        % Store network Effect of EB Plans:
        AE_exp_net_orig = AE_exp_orig(:,3);
        
        % Determine length of MonteCarlo:
        nMC = 1e5;
        
        % Store results into a benchmark array:
        original = [ mean(AE_tax_net_orig)/std(AE_tax_net_orig)  ...
            1-sum(AE_tax_net_orig<0)/nMC; ...
            mean(AE_exp_net_orig)/std(AE_exp_net_orig) ...
            1-sum(AE_exp_net_orig<0)/nMC];
        
        % Remove posterior from the workspace
        clearvars posterior
        
        
        %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        % Carry out a Placebo Analysis: simulate the network by resambling
        % the original one and re-run the placebo estimation analysis by
        % calculating the network effect of TB and EB plans:
        
        % Make sure you are running a baseline analysis:
        if strcmp(inversion,"yes")
            error("Placebo can't run with inverted model")
        elseif strcmp(model,"dynamic")
            error("Placebo can't run with dynamic model")
        end
        
        
                
        % Initialize: (first floor TB, second floor is EB)
        placebo_results = zeros(nPl,2,2);
        
        % Run Placebo
        parfor pl = 1 : nPl
            
            
            %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            % Simulate the network:
            [Adown_syn,~,ub_down_syn] = network_simulation(Adown);
            [Aup_syn,~,ub_up_syn] = network_simulation(Aup);
            % Store data into cell array:
            data_MCMC_pl = {dy_m dy ...
                d_down d_up t_down t_up ...
                k ...
                Adown_syn Aup_syn 0 0 ub_down_syn ub_up_syn ...
                X X_down X_up};
            
            
            
            
            %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            % Estimate parameters via Bayesian MCMC:
            posterior = Bayesian_MCMC(data_MCMC_pl,nMCMC,n,T);
            
            
            
            
            %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            % Calculate Average Effects:
            % Simulation of average effects of a TB plan:
            AE_tax = ave_eff_MC(posterior.rho_down,...
                posterior.beta(:,index_tau),n,Adown_syn,"tax",data,uncertainty,...
                share_TB,industry_weight_string,data.ind_weight_TB,periods);
            % Store network Effect of TB Plans:
            AE_tax_net = AE_tax(:,3);
            
            
            
            % Simulation of average effects of an EB plan:
            AE_exp = ave_eff_MC(posterior.rho_up,...
                posterior.beta(:,index_gamma),n,Aup_syn,"exp",data,uncertainty,...
                share_EB,industry_weight_string,data.ind_weight_EB,periods);
            % Store network Effect of EB Plans:
            AE_exp_net = AE_exp(:,3);
            
            % Calculate Descriptive statistics we are interested in, and
            % store the results into a 3D array: 1) placebo results ; 2)
            % there are two columns, the first one is the mean/std of the
            % posterior of the network effect, the second one is the
            % Pr(X<0); 3) first floor is TB, second floor is EB.
            placebo_results(pl,:,:) = cat(3, ...
                [mean(AE_tax_net)/std(AE_tax_net) 1-sum(AE_tax_net<0)/nMC], ...
                [mean(AE_exp_net)/std(AE_exp_net) 1-sum(AE_exp_net<0)/nMC]);
            
        end
        
        
        % Plot Results of Placebo Analysis:
        fig_pl = figure_placebo(placebo_results,original,periods);
        
        
end



%% Industry Specific Network Effects:

if strcmp(model,"static") && strcmp(placebo,"no")
    
    
    % Construct Customerness (or "Downstreamness"):
    Customerness = sum( inv(eye(n)-Adown)-eye(n) ,2);
    
    % Construct Supplierness (or "Upstreamness"):
    Supplierness = sum( inv(eye(n)-Aup)  -eye(n) ,2);
    
    % Construct Average Network Industry Specific Response:
    mean_NEi_TB = mean(posterior.AE_tax_ind(:,:,3),2);
    mean_NEi_EB = mean(posterior.AE_exp_ind(:,:,3),2);
    
    % Construct a figure which plots: 1) the mean, industry specific 
    % network effect of TB plans against the "Customerness"; 2) the mean,
    % industry specific network effect of EB plans against the
    % "Supplierness":
    data_industry_analysis = { Customerness mean_NEi_TB ...
        Supplierness mean_NEi_EB data.industries};
    
    % Plot:
    fig_industry_effects(data_industry_analysis,"TB");
    
    
    
    
end



%% Plot Fiscal Plans:

[fig1,fig2] = figure_fiscal_plans(data,tb,eb);



%% Partitioning Table:

table_partitioning = partitioning(A,A_hat_tr);


%% Plot Kernel of a Non-Standardized Beta Distribution:

% Support:
a = 0;
b = ub_down;

% Shape parameter:
d = 1.1;

% Plot:
fig_beta_kernel = gen_beta_plot(a,b,d);



%% Export Tables to Excel:

% Export the "core parameters estimates" Table to Excel:
writetable(table_IIIa(1:10,:),'tables.xlsx','Sheet','Parameters')

% Export theAverage Effects Table to Excel:
writetable(tableIIIb,'tables.xlsx','Sheet','Average Effects')




%% Plot Acceptance Rate:


subplot(2,1,1),plot(posterior.acc_rate(:,1),'b','Linewidth',1.5)
grid on
subplot(2,1,2),plot(posterior.acc_rate(:,2),'b','Linewidth',1.5)
grid on



