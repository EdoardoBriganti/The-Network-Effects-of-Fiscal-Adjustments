function [fisc_cons_ols,fisc_cons_summary,fiscal_data,var] = fiscal_shock_analysis(data)
%{

===========================================================================
PURPOSE: 
Provides summary statistics of the fiscal shocks. 

===========================================================================
 INPUT   'data'    is the structure which contains all the relevant time
                   series in our dataset. 'data' is uploaded through the
                   upload function.
            'T'    is the length of the tme series (T=37)  

OUTPUT  'summary'   is a  5x5x2 matrix, which contains the summary of the
                    fical shocks. The first "floor" contains data referring
                    to taxes, the second floor contains the data for the
                    spending shocks. Rows number: 1) number of non zero
                    entries; 2) mean (of non zero entries); 3) std (of non
                    zero entries); 4) min; 5) max. The columns: 1)
                    unanticiapted; 2) anticipated; 3) future one year
                    ahead; 4) future, two years ahead; 5) future, 3 years
                    ahead.
            'b'     is a 2x3x2 array which contains the result of the
                    regressions. The first "floor" contains data referring
                    to taxes, the second floor contains the data for the
                    spending shocks. First row is constant, second row is
                    the regression coefficient. First column is regression
                    of SA1 on SU, and so on.
   'fiscal_data'    It is a table which represents the fiscal shock
                    database that we employ in our analysis. 

%=========================================================================%

Favero, Karamysheva, Briganti 
"The network effects of fiscal adjustments" (2019)

%=========================================================================%
%}

x = cat(3, [data.SUT data.SAT data.SAT1 data.SAT2 data.SAT3 ],...
    [data.SUG data.SAG data.SAG1 data.SAG2 data.SAG3]);
k = size(x,2);
fisc_cons_summary = zeros(5,k,2);
fisc_cons_ols = zeros(2,k,2);
for type = 1:2
    
    
    X = x(:,1,type);
    % Remove zeros;
    ind = find(X);
    X(X==0) = [];
    X = [ones(length(X),1) X];
    for time = 1:k
        
        s = x(:,time,type);
        fisc_cons_summary(1,time,type) = sum(s~=0);
        s(s==0) =[ ];
        fisc_cons_summary(2,time,type) = mean(s);
        fisc_cons_summary(3,time,type) = std(s);
        fisc_cons_summary(4,time,type) = min(s);
        fisc_cons_summary(5,time,type) = max(s);
        % Regress shocks on unanticipated components conditional on those
        % years where an unanticipated shock occurs.
        Y = x(:,time,type);
        fisc_cons_ols(:,time,type) = X\Y(ind,:);
        % Simulate the variance of the regression:
        t = size(X,1);
        e = Y(ind,:) - X * fisc_cons_ols(:,time,type);
        var(time,type) = ( e' * e ) / (t - 2);
    end
end

fiscal_data = x;

end