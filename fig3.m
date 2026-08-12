[c, ch] = radial_coupled_pde_model(1/10, 10, 10, 0.25, 0.99*0.75, 250, 1e-6, logspace(-3, 1, 5));

% Plot concentrations
for i = 1 : 5
	plot(linspace(0, 1, 251), c(i, :));
	hold on;
end

for i = 1 : 5
	plot(linspace(0, 1, 251), ch(i, :));
	hold on;
end

% Plot contour concentrations
radii = linspace(0, 1, 251);

X = zeros(251, 251);
Y = X;

for i = 1 : 201
    X(i, :) = radii(i)*cos(linspace(0, 2*pi, 201));
    Y(i, :) = radii(i)*sin(linspace(0, 2*pi, 201));
end

for i = 1 : 5
	contourf(X, Y, repmat(c(i, :)', [1 251]), 200, 'edgecolor', 'none');
	hold on;
end

% Plot R comparisons
for R = [0.1 0.5 1 5 10 50 100 1000]
    [c, ch] = radial_coupled_pde_model(1/10, 10, R, 0.25, 0.99*0.75, 250, 1e-6, logspace(-3, 1, 200));
    
    semilogx(logspace(-3, 1, 200), c(:, 125));
    hold on;
end
