function table = table_average_effect(posterior)

% Posterior Distributions of Average Effects:
% TB Plans:
ATE_TB = posterior.AE_tax(:,1);
AIE_TB = posterior.AE_tax(:,2);
ANE_TB = posterior.AE_tax(:,3);
% Eb Plans:
ATE_EB = posterior.AE_exp(:,1);
AIE_EB = posterior.AE_exp(:,2);
ANE_EB = posterior.AE_exp(:,3);

% Store Average Effects in a matrix:
data = [ATE_TB AIE_TB ANE_TB ATE_EB AIE_EB ANE_EB];

% Calculate Length of the Distribution (Number of observations)
MC  = size(data,1);

% Preallocate Table:
tab = zeros(10,size(data,2));

% Construct Descriptive statistics 
for j = 1 : size(data,2)
    
    % Expected Value:
    tab(1,j) = mean(data(:,j));
    % Standard Deviation:
    tab(2,j) = std(data(:,j));
    % Probability of being negative:
    tab(3,j) = sum(data(:,j)<0)/MC;
    
    % Percentiles:
    data(:,j) = sort(data(:,j));
    tab(4,j) = data(round(0.05*MC),j);
    tab(5,j) = data(round(0.1*MC),j);
    tab(6,j) = data(round(0.16*MC),j);
    tab(7,j) = data(round(0.5*MC),j);
    tab(8,j) = data(round(0.84*MC),j);
    tab(9,j) = data(round(0.9*MC),j);
    tab(10,j) = data(round(0.95*MC),j);

end

% Calculate Share of Total Effect coming from Instantaneous and Network
% Effects, on average (calculated on the expected value):
div   = [repmat(tab(1,1),[1,3]) repmat(tab(1,4),[1,3])];
tab   = [tab(1,:) ; tab(1,:) ./ div ; tab(2:end,:)];

% Construct Table:
table = array2table(tab,'VariableNames',{'TB_tot','TB_instant',...
    'TB_network','EB_tot','EB_instant','EB_network'});


end