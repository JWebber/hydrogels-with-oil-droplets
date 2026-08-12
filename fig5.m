K = 10;

for psi = [0 0.05 0.1 0.25 0.5]
    [c, f] = coupled_pde_model(1, K, 100, psi, 0.99*(1-psi), 100, 2e-5, logspace(-2, log10(50), 100));
    semilogx(logspace(-2, log10(50), 100), f*(1-psi));
    hold on;
end

for psi = [0 0.05 0.1 0.25 0.5]
    [c, f] = effective_model(1, 0.5, K, 100, psi, 0.99*(1-psi), 100, 2e-5, logspace(-2, log10(50), 100));
    semilogx(logspace(-2, log10(50), 100), f, 'linestyle', '--');
    hold on;
end