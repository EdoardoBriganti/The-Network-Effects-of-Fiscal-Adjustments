function fig_industry_effects = fig_industry_effects(data_industry_analysis,...
    industry_effects2plot)

% Extract Data:
Customerness = cell2mat(data_industry_analysis(1));
mean_NEi_TB  = cell2mat(data_industry_analysis(2));
Supplierness = cell2mat(data_industry_analysis(3));
mean_NEi_EB  = cell2mat(data_industry_analysis(4));
industries   = data_industry_analysis{5};

% Order of original data:
% data_industry_analysis = { Customerness   mean_NEi_TB ...
%         Supplierness    mean_NEi_EB    data.industries};


fig_industry_effects = figure;
switch industry_effects2plot
    case "both"
        
        % TB:
        subplot(1,2,1),scatter(Customerness,mean_NEi_TB,17,'b')
        lsline
        text(Customerness+0.01,mean_NEi_TB+0.01,industries)
        grid on
        xlabel("Customerness")
        ylabel('$\mathbb{E}(\textit{Network Effect}_i$ of TB Plans','Interpreter','Latex')
        set(gca,'Fontsize',25)
        % EB:
        subplot(1,2,2),scatter(Supplierness,mean_NEi_EB,17,'b')
        lsline
        text(Supplierness+0.01,mean_NEi_EB+0.01,industries)
        grid on
        xlabel("Supplierness",'Interpreter','Latex')
        ylabel('$\mathbb{E}(\textit{Network Effect}_i$ of EB Plans','Interpreter','Latex')
        set(gca,'Fontsize',25)
        
    case "TB"
        % TB:
        scatter(Customerness,mean_NEi_TB,17,'b')
        lsline
        text(Customerness+0.01,mean_NEi_TB+0.01,industries)
        grid on
        xlabel("Customerness",'Interpreter','Latex')
        ylabel('$\bf{E}( Network Effect )_i$ of TB Plans','Interpreter','Latex')
        set(gca,'Fontsize',25)
        
    case "EB"
        % EB:
        scatter(Supplierness,mean_NEi_EB,17,'b')
        lsline
        text(Supplierness+0.01,mean_NEi_EB+0.01,industries)
        grid on
        xlabel("Supplierness",'Interpreter','Latex')
        ylabel('$\mathbb{E}$ of EB Plans','Interpreter','Latex')
        set(gca,'Fontsize',25)
end




end