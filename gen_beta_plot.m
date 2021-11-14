function fig = gen_beta_plot(a,b,d)


% Shape and Scale parameter:
A = d;
B = d;

% Vector of X:
x = transpose(linspace(a,b,1e3));

% Density of Generalize (or non std.) Beta Distribution:
f = ( (x-a).^(A-1) .* (b - x).^(B-1) ) ./ ( beta(A,B) * (b-a)^(A+B-1) );


fig = figure;
plot(x,f,'Linewidth',4)
grid on
xlabel('$\rho$','Interpreter','Latex')
ylabel('Density','Interpreter','Latex')
text(b+0.01,a-0.015,'$\lambda_{max}^{-1}(A)$','Interpreter','Latex','Color','red','Fontsize',25)
set(gca,'Fontsize',25)

end



