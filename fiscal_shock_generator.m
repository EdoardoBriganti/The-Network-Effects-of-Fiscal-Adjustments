function [tax,exp] = fiscal_shock_generator(fisc_cons_ols,fisc_cons_summary,uncertainty)
%{

===========================================================================
PURPOSE:
This function construct the fiscal shocks, that is, given a 1% of GDP
fiscal shock, how it is distributed among its components (unanticipated,
and futures)? 
There are to main ways to do that: 
1) uncertainty: 'no'. Under this option, the function distributes the shock
                      using the mean value of the time series.
2) uncertainty: 'yes'. Under this option, the function randomly draw the
                       values of the unanticipated shock from a uniform 
                       distribution on the interval constructed from the
                       mean value plus/minus a std. 
!!) The shocks are normalized to one. For the purposes of our analysis we
    don't give a shock to the anticipated component. For this reason the
    second entry of 'tax' and 'exp' are zeros.

===========================================================================
INPUT:
          b: is a 3D-array of regression coefficients. We regress the
             unexpected components over the remaining components. For
             further reference, this array is the output of the function
             called "fiscal_shock_analysis".
    summary: is a matrix of descriptive statistics of the fiscal shocks
             database. For further reference, this array is the output of 
             the function called "fiscal_shock_analysis".

 fiscal_shock_uncertainty: this is a string variable. It can be either
                           'yes' or 'no'. If equal to 'yes', the function
                           generates shocks from a distribution (see
                           below). If equal to 'no' the function uses a
                           shock distribution which resembles the orginal
                           database (average shock conditional on shock
                           occurrence).

OUTPUT:  
        tax: is 5x1 vector which contains the shock distribution normalized
             to 1.
        exp: is 5x1 vector which contains the shock distribution normalized
             to 1. In this case, the unexpected component, when drawn from
             the distribution (because the fiscal_shocks_uncertainty is
             'yes') could be negative, if that happen, the function
             automatically sets it equal to zero.

===========================================================================


%}



% Step 1) define the parameters of the uniform distribution from which we
% draw the unexpected components, for both taxes and expenditures.
mu_tax    = fisc_cons_summary(2,1,1);
std_tax   = fisc_cons_summary(3,1,1);
mu_exp    = fisc_cons_summary(3,1,2);
std_exp   = fisc_cons_summary(2,2,2);

% Step 2)
% Draw unexpected components from the uniform distribution whose support is
% an interval of 2 std range centerd on the mean (both mean and std are
% calculated conditional on those years where an unexpected shock occur).
% If it is required to calculate the shocks without uncertainty, then the
% mean values of every variable is considered for the calculation of the
% distribution of the simulated shock.
S = zeros(5,2);
if strcmp(uncertainty,'yes')==1
    S(1,1)    = (mu_tax - std_tax) + (mu_tax + std_tax) * rand;
    S(1,2)    = (mu_exp - std_exp) + (mu_exp + std_exp) * rand;
else
    S(1,1) = mu_tax;
    S(1,2) = mu_exp;
end

% Step 3)
% After drawing the unexpected component, use the regression analysis
% performed by the function "shocks" to calculate the corresponding future
% components. 
for type = 1:2 % tax or exp?
    for time = 1:3  % what kind of future adjustment
        S(time+2,type) = [1 S(1,type)] * fisc_cons_ols(:,time+2,type);
    end
end

% Step 4)
% Normalize the shocks, so that they sum up to one. 
S = S ./ sum(S);

% Step 5)
% Store results into two OUTPUT vectors and deal with potential negative
% shocks: 
% Step 5a) if the unexpected component is negative aafter normalization,
% replace the unexpected component with zero and use the constant as shocks
% for the remaining components. 
tax = S(:,1);
if tax(1) < 0
    tax(1:2,1) = zeros(2,1);
    tax(3:end,1) = squeeze(fisc_cons_ols(1,3:end,1))';
end
exp = S(:,2);
if exp(1) < 0
    exp(1:2,1) = zeros(2,1);
    exp(3:end,1) = squeeze(fisc_cons_ols(1,3:end,2))';
end
% Step 5b) After having rule out negative unexpected component, we still
% need to rule out potential negative values for the future components. If
% one of those values are negative, replace them with zero, and then
% perform a normalization.
tax(tax<0) = 0;
tax = tax ./ sum(tax);
exp(exp<0) = 0;
exp = exp ./ sum(exp);

%===================================================%
% If you want to use the distribution of the shocks that we used to employ
% in the former versions of the paper, set the string "uncertainty" to be
% equal to 'old':

if strcmp(uncertainty,'old')==1
    tax = [.115 ; 0 ; .885 ; 0 ; 0];
    exp = [.198 ; 0 ; .802  ; 0 ; 0];
end


end