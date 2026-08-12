[X, Y] = meshgrid(linspace(0, 1, 250), logspace(-3, 3, 250));

contourf(X, Y, log10(1./(1-Y.*X./(0.99.*(1-X)+Y.*X))), 250, 'edgecolor', 'none');
set(gca, 'yscale', 'log');
clim([0 2]);