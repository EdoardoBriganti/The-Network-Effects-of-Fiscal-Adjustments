function fig = cumIRF(posterior)

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Step 1: Construct Average Cumulative Impulse Responses:

% TB Plans:
cum_irf_tax_tot = squeeze(mean(posterior.cum_irf_tax_tot,1));
cum_irf_tax_dir = squeeze(mean(posterior.cum_irf_tax_dir,1));
cum_irf_tax_net = squeeze(mean(posterior.cum_irf_tax_net,1));
% EB Plans:
cum_irf_exp_tot = squeeze(mean(posterior.cum_irf_exp_tot,1));
cum_irf_exp_dir = squeeze(mean(posterior.cum_irf_exp_dir,1));
cum_irf_exp_net = squeeze(mean(posterior.cum_irf_exp_net,1));



%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Step 2: Concatenate CumIRF into a 3D matrix:

% TB Plans:
cum_irf_tax = cat(3,cum_irf_tax_tot,cum_irf_tax_dir,cum_irf_tax_net);
% EB Plans:
cum_irf_exp = cat(3,cum_irf_exp_tot,cum_irf_exp_dir,cum_irf_exp_net);

% !!!) Order: 1) Total ; 2) Direct; 3) Network



%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Step 3: Calculate Median and Std of each CumIRF of each MC distribution:

% Construct Indices of percentiles:
MC = size(cum_irf_tax,2);
p05 = ceil(0.05*MC);
p50 = ceil(0.5*MC);
p95 = ceil(0.95*MC);


% Select Objects to plot:
% TB Plan:
cum_irf_tax = sort(cum_irf_tax,2); 
cum_irf_tax_plot = cat(2,cum_irf_tax(:,p05,:),cum_irf_tax(:,p50,:),cum_irf_tax(:,p95,:));


% EB Plan:
cum_irf_exp = sort(cum_irf_exp,2); 
cum_irf_exp_plot = cat(2,cum_irf_exp(:,p05,:),cum_irf_exp(:,p50,:),cum_irf_exp(:,p95,:));

y_lim = min([squeeze(min(min(cum_irf_tax_plot,[],1),[],2)) ...
    squeeze(min(min(cum_irf_exp_plot,[],1),[],2))],[],2);



%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% Step 4: Plot Cumulative Impluse Response with 90% Confidence Bands:

str_title = ["Total" ; "Instantaneous" ; "Network"];
horz = size(cum_irf_tax,1)-1;
X    = transpose(0:1:horz);
fig  = figure(1);
for i = 1:3
    
    % TB Plan:
    subplot(2,3,i),shade(X,cum_irf_tax_plot(:,1,i),X,cum_irf_tax_plot(:,3,i),'FillType',[1 2;2 1]);
    hold on
    subplot(2,3,i),plot(X,cum_irf_tax_plot(:,1,i),'--r',...
        X,cum_irf_tax_plot(:,2,i),'b-d',X,cum_irf_tax_plot(:,3,i),'--r','Linewidth',2)
    grid on
    xlabel('Years','Interpreter','Latex')
    ylabel('\% of GDP','Interpreter','Latex')
    title("TB Plan - " + str_title(i),'Interpreter','Latex')
    set(gca,'Fontsize',24)
    %ylim([y_lim(i) , 0])
    hold off
    
    % EB Plan
    subplot(2,3,i+3),shade(X,cum_irf_exp_plot(:,1,i),X,cum_irf_exp_plot(:,3,i),'FillType',[1 2;2 1]);
    hold on
    subplot(2,3,i+3),plot(X,cum_irf_exp_plot(:,1,i),'--r',...
        X,cum_irf_exp_plot(:,2,i),'b-d',X,cum_irf_exp_plot(:,3,i),'--r','Linewidth',2)
    grid on
    xlabel('Years','Interpreter','Latex')
    ylabel('\% of GDP','Interpreter','Latex')
    title("EB Plan - " + str_title(i),'Interpreter','Latex')
    set(gca,'Fontsize',24)
    %ylim([y_lim(i) , 0])
    hold off
end


end