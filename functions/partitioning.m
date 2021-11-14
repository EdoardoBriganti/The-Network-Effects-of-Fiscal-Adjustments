function tab = partitioning(A,A_hat_tr)

% Size of the network:
n = size(A,1);
% Assume no prior knowledge on the estimate of rho:
rho = 1;

% Calculate Lenotief Inverse (total feedback effect)
H        = inv( eye(n) - rho*A);
H_hat_tr = inv( eye(n) - rho*A_hat_tr);

% Partition the effect:
max_order = 5;
downstream = zeros(max_order+1,1);
upstream   = zeros(max_order+1,1);

for j = 0 : max_order
    
    % Downstream
    downstream(j+1,1) = mean( ((rho*A)^j * ones(n,1) ) ./ (H * ones(n,1) ) );
    
    % Upstream
    upstream(j+1,1)   = mean( ((rho*A_hat_tr)^j * ones(n,1) ) ./ (H_hat_tr * ones(n,1) ) );
    
    % String of Effects:
    order(j+1,1) = "Order " + string(j);
     
end
% Calculate cumulative effect:
downstream_cum = cumsum(downstream);
upstream_cum   = cumsum(upstream);
% Construct Table:
tab = table(order,downstream,downstream_cum,upstream,upstream_cum);



end