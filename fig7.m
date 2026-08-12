K = 1;

for S = [0.01 0.1 1 10 100]
    [c, ch, f] = spherical_coupled_pde_model(1, K, S, 0.25, 0.99*(1-0.25), 100, 2e-5, logspace(-2, 2, 100));

    r = linspace(0, 1, 100+1);
    m_r = repmat(0.5*(r(2:end)+r(1:end-1)), [100, 1]);

    m_c = 0.5*(c(:, 2:end)+c(:, 1:end-1));
    m_cH = 0.5*(ch(:, 2:end)+ch(:, 1:end-1));

    semilogx(logspace(-2, 2, 100), sum(4*pi*m_r.^2.*m_c, 2)/100 + sum(4*pi*m_r.^2.*m_cH, 2)/100);
end

[c, f] = spherical_effective_model(1, K, 1, 0.25, 0.99*(1-0.25), 100, 2e-5, logspace(-2, 2, 100));

r = linspace(0, 1, 100+1);
m_r = repmat(0.5*(r(2:end)+r(1:end-1)), [100, 1]);

m_c = 0.5*(c(:, 2:end)+c(:, 1:end-1));

semilogx(logspace(-2, 2, 100), 2*sum(4*pi*m_r.^2.*m_c, 2)/100, 'k--');