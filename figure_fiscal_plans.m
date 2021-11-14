function [fig1,fig2] = figure_fiscal_plans(data,tb,eb)

% Figure 1: Fiscal Adjustment Database:
data2plot1 = [ data.SU.*tb data.SA.*tb data.SA_F.*tb ...
    data.SU.*eb data.SA.*eb data.SA_F.*eb  ];

J = size(data2plot1,2);

tit = ["$f_t^u \cdot TB_t$" "$f_t^a \cdot TB_t$" "$f_t^f \cdot TB_t$" ...
    "$f_t^u \cdot EB_t$" "$f_t^a \cdot EB_t$" "$f_t^f \cdot EB_t$"  ];

fig1 = figure;
for j = 1:J
    
   subplot(2,3,j),plot(data.years,data2plot1(:,j),'b','Linewidth',2)
   xlabel('Years','Interpreter','Latex')
   ylabel('\% of GDP','Interpreter','Latex')
   title(tit(j),'Interpreter','Latex')
   grid on
   set(gca,'Fontsize',25)
   
end


% Figure 2: labelling
tax = sum(data.tax_shocks,2);
tax(tax<0) = 0;
exp = sum(data.exp_shocks,2);
f_tax = tax ./ (tax+exp);
f_exp = exp ./ (tax+exp);

fig2 = figure;
subplot(2,1,1),bar(data.years,[f_tax f_exp].*tb,'Stacked')
grid on
ylabel('\%','Interpreter','Latex')
title('TB Fiscal Adjustment composition','Interpreter','Latex')
legend("Tax Rise","Spending Cut",'Location','NorthEast')
set(gca,'Fontsize',25)

subplot(2,1,2),bar(data.years,[f_tax f_exp].*eb,'Stacked')
grid on
ylabel('\%','Interpreter','Latex')
title('EB Fiscal Adjustment composition','Interpreter','Latex')
legend("Tax Rise","Spending Cut",'Location','NorthWest')
set(gca,'Fontsize',25)


end