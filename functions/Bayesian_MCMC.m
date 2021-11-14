function posterior = Bayesian_MCMC(data_MCMC,nMCMC,n,T)




%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Extract Data:
dy_m   = cell2mat(data_MCMC(1));

d_down = cell2mat(data_MCMC(2));
d_up   = cell2mat(data_MCMC(3));
t_down = cell2mat(data_MCMC(4));
t_up   = cell2mat(data_MCMC(5));

k      = cell2mat(data_MCMC(6));

Adown  = cell2mat(data_MCMC(7));
Aup    = cell2mat(data_MCMC(8));

ub_up   = cell2mat(data_MCMC(9));
ub_down = cell2mat(data_MCMC(10));

X      = cell2mat(data_MCMC(11));

% Necessary elements to implement MH for rho1:
X_down = cell2mat(data_MCMC(12));
dy_down = dy_m(:,find(d_down')); % (n x t_down) matrix

% Necessary elements to implement MH for rho2:
dy_up  = dy_m(:,find(d_up')); % (n x t_up) matrix
X_up   = cell2mat(data_MCMC(13));

% This is how data were saved: order must be consistent.
% data_MCMC = {dy_m  ...
%     d_down  d_up   t_down  t_up ...
%     k ...
%     Adown  Aup   ub_up   ub_down  ...
%     X  X_down  X_up};



%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Set up prios:

% 1) beta~MVN(c,L)  ====> beta|Data,Teta~MVN(c*,L*)
prior.c = zeros(k,1);
prior.L = diag(1e12*ones(k,1));     % Diffuse Prior on the beta parameter (k x k) matrix.
Inv_L   = diag( 1./diag(prior.L) ); % Invert Diagonal Matrix

% 2) sig2~IG(a,b) ====> sig2|data,Teta~IG(a*/2,b*/2)
%    a* = nT+2a
%    b* = M + 2b
%    if a = b = 0 ====> we are putting a Jefferey's prior on sig2.
prior.a = 0;
prior.b = 0; 

% 3) r/vi~Chi2(r) ====> vi|Data,Teta~IG(q1/2,q2/2)
%    q1 = r+T
%    q2 = (1/sig2)*F{i} + r
prior.r = 3;

% 4) Draw using Metropolis Hastings the values for rho1 and rho2.
%    we put the same beta prior on both rho1 and rho2, to resamble a 
%    uniform distribution with less density on extreme values of rho:
%    rho1~beta(d,d) rho2~beta(d,d)
prior.d = 1.1;





%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Set up initial values:
c_prop1  = 0.1; % Variance of the proposal distribution of rho1.
c_prop2  = 0.1; % Variance of the proposal distribution of rho2.
vi       = ones(n,1);
V        = diag(vi);
sig2     = 1;
rho_down = 0.5;
rho_up   = 0.5;
beta     = zeros(k,1);

% These initial values will be updated step by step within the MH algorithm:
Inv_V = diag( 1./diag(V) );

% Don't vectorize to avoid inversion of large matrices which slow down the
% code:
for t = 1:T
    
    % Networked-Purged Dependent Variable:
    Z(:,t) = dy_m(:,t) - rho_down * (Adown * dy_m(:,t)) * d_down(t) - ...
        rho_up *  (Aup * dy_m(:,t)) * d_up(t) ; % (n x T) matrix
    
    % Residuals:
    eps(:,t) = Z(:,t) - X(:,:,t)*beta;  % (n x T) matrix
    
    % Sum of Squared Residuals addend:
    M_m(t) = eps(:,t)' * Inv_V * eps(:,t);  % (T x 1) vector
    
    % Sum of Squared Errors by unit:
    F_m(:,t) = eps(:,t).^2;  % (n x T) matrix
    
end

% Sum of Squared Residuals:
SSR = sum(M_m);

% Sum of Squared Errors by unit:
eps2 = sum(F_m,2);




%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% STEP 3) Bayesian MCMC:

% Simulation length
niter = floor(nMCMC/0.9);

% Preallocate:
post_beta = zeros(k,niter);
post_sig2 = zeros(niter,1);
post_vi   = zeros(n,niter);
post_rho_down = zeros(niter,1);
post_rho_up = zeros(niter,1);
acceptance = zeros(niter,2);
acceptance_rate_rho_down = 0;
acceptance_rate_rho_up = 0;

% Start the Bayesian Markov Chain Monte Carlo Simulation
for mcmc = 1:niter
    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % STEP 1) Gibbs Sampling:
    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
    % 1a) - Draw beta:
    
    % Construct useful elemnts:
    xVx = zeros(k,k);
    xVz = zeros(k,1);
    for t = 1:T
        % Exploit the fact that V is diagonal to save computational power:
        %xV  = X(:,:,t)' ./ vi';
        % Calculate useful elements:
        xVx = X(:,:,t)' * Inv_V * X(:,:,t) + xVx;
        xVz = X(:,:,t)' * Inv_V * Z(:,t) + xVz;
    end
    
    % Variance of Conditional Full Likelihood of beta
    L_star = sig2 * inv( xVx + sig2*Inv_L );
    
    % Mean of Conditional Full Likelihood of beta
    c_star = ( xVx + sig2*Inv_L ) \ ( xVz + sig2*Inv_L*prior.c );
    
    % Draw from MV-Normal N=
    post_beta(:,mcmc) = mvnrnd(c_star,L_star);
    
    % Update beta
    beta   = post_beta(:,mcmc);  
    
    % Update elements which depend on beta:
    SSR = 0;
    eps2 = zeros(n,1);
    for t = 1:T
        % Residuals:
        eps(:,t) = Z(:,t) - X(:,:,t)*beta;  % (n x T) matrix
    
        % Sum of Squared Residuals addend:
        eVe = eps(:,t)' * Inv_V * eps(:,t);  % (T x 1) vector
        SSR = SSR + eVe;
        
        % Sum of Squared Errors by unit:
        F_m(:,t) = eps(:,t).^2;  % (n x T) matrix
        
        % Squared of errors by unit:
        eps2 = eps2 + eps(:,t).^2;
    end
    
    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
    % 1b) Draw sigma:
    
    % Parameter of the inverse gamma:
    a_star = n*T + 2*prior.a;
    b_star = SSR + 2*prior.b;
    
    % Draw Sigma2:
    post_sig2(mcmc,1) = b_star/chi2rnd(a_star);
    % Update Sigma2
    sig2 = post_sig2(mcmc,1);  % scalar
    
    
    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=%
    % 1c) Draw vi:
    
    % Parameters of the Inverse Gamma:
    q1 = prior.r + T;
    q2 = (1/sig2).*eps2 + prior.r;
    
    % Draw vi for each industry:
    post_vi(:,mcmc) = q2./chi2rnd(q1,[n,1]);
    % Update vi
    vi = post_vi(:,mcmc);     % (n x 1) vector  

    % Construct inverse of V:
    Inv_V = diag(1./vi);    % (n x n) matrix

    
    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    % STEP 2) Metropolis-Hastings:
    
    % Step 2a - Draw rho1:
    [rho_down,c_prop1,acceptance_rate_rho_down] = metropolis_hastings(rho_down,...
        c_prop1,acceptance_rate_rho_down,n,Adown,t_down,sig2,dy_down,Inv_V,beta,X_down,prior.d,mcmc,0,ub_down);
    post_rho_down(mcmc,1) = rho_down;
    
    % Step 2b - Draw rho2:
    [rho_up,c_prop2,acceptance_rate_rho_up] = metropolis_hastings(rho_up,...
        c_prop2,acceptance_rate_rho_up,n,Aup,t_up,sig2,dy_up,Inv_V,beta,X_up,prior.d,mcmc,0,ub_up);
    post_rho_up(mcmc,1) = rho_up;
    
    % Update elements which depend on rho1 and rho2: Z, the Networked-Purged Dependent Variable
    for t = 1:T
        Z(:,t) = dy_m(:,t) - rho_down * (Adown * dy_m(:,t)) * d_down(t) - ...
            rho_up *  (Aup * dy_m(:,t)) * d_up(t) ; % (n x T) matrix
    end
    
    acceptance(mcmc,:) = [acceptance_rate_rho_down acceptance_rate_rho_up];
end

% Burn-in (remove the first 10% of simulated data)
burn_in = floor(.1 * niter);

% Store results
posterior.rho_down = post_rho_down(burn_in:end,:);
posterior.rho_up   = post_rho_up(burn_in:end,:);
posterior.beta     = transpose(post_beta(:,burn_in:end));
posterior.sig2     = post_sig2(burn_in:end,:);
posterior.vi       = post_vi(:,burn_in:end);
posterior.acc_rate = acceptance;

end