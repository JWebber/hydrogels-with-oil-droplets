Ks = linspace(0.01, 100, 100);
load('fit_diffs.mat');

fits = real(fits);

fitfun = @(p)((p(1) + p(2).*Ks).^(-1) - fits');

fit_params = lsqnonlin(fitfun, [1, 1]);

plot(Ks, fits, 'r.');
hold on;
plot(Ks, fitfun(fit_params)+fits');

matlab2tikz('fitted_diffs.tikzdata');