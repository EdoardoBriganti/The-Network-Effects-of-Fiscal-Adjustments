function fig_pl = figure_placebo(placebo_results,original,periods,placebo_type)


switch periods
    case "1year"
        
        % Plot only the network effect:
        % TB Plans:
        fig_pl = figure; 
        subplot(1,2,1),scatter(placebo_results(:,1,1),placebo_results(:,2,1),17,'b')
        hold on
        subplot(1,2,1),scatter(original(1,1),original(1,2),80,'r','MarkerFaceColor','r')
        xlabel('Standardized Mean of Average Network Effect','Interpreter','Latex')
        ylabel('$Pr(X<0)$','Interpreter','Latex')
        legend('Simulated Spatial Matrix','Original Spatial Matrix','location','northwest')
        grid on
        title('Placebo Test - TB plan','Interpreter','Latex')
        set(gca,'Fontsize',25)
        hold off
        % EB Plans:
        subplot(1,2,2),scatter(placebo_results(:,1,2),placebo_results(:,2,2),17,'b')
        hold on
        subplot(1,2,2),scatter(original(2,1),original(2,2),80,'r','MarkerFaceColor','r')
        xlabel('Standardized Mean of Average Network Effect','Interpreter','Latex')
        ylabel('$Pr(X<0)$','Interpreter','Latex')
        legend('Simulated Spatial Matrix','Original Spatial Matrix','location','northwest')
        grid on
        title('Placebo Test - EB plan','Interpreter','Latex')
        set(gca,'Fontsize',25)
        
    case "2years"
        
        % Plot only the network effect of TB plans:
        fig_pl = figure; 
        scatter(placebo_results(:,1,1),placebo_results(:,2,1),80,'b','MarkerFaceColor','y')
        hold on
        scatter(original(1,1),original(1,2),150,'k','MarkerFaceColor','r')
        xlabel('Standardized Mean of $\textit{ANE}_{TB}$','Interpreter','Latex')
        ylabel('$Pr(\textit{ANE}_{TB}>0)$','Interpreter','Latex')
        legend('Simulated Spatial Matrix','Original Spatial Matrix','location','northwest')
        grid on
        title(placebo_type,'Interpreter','Latex')
        set(gca,'Fontsize',30)
        hold off
end




end