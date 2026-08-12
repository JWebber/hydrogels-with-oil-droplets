function [C, CH, F] = spherical_coupled_pde_model(D, K, RS, Psi, Vw, Nx, dt, times)
    sz = max(size(times, 1), size(times, 2));

    t = 0;
    reportIndex = 1;

    C = zeros(sz, Nx+1);
    CH = C;
    F = zeros(sz, 1);

    % BUILD UP A NONTRIVIAL INITIAL CONDITION
    c = ones(1, Nx+1);
    cHat = K*ones(1, Nx+1);

    r = linspace(0, 1, Nx+1);
    m_r = 0.5*(r(2:end)+r(1:end-1));

    while (reportIndex <= sz)
        if(t >= times(reportIndex))
            CH(reportIndex, :) = cHat;
            C(reportIndex, :) = c;
            F(reportIndex) = -D*(1-Psi)*Nx*(c(end)-c(end-1));
            reportIndex = reportIndex + 1;
        end

        m_d_c = Nx*(c(2:end)-c(1:end-1)).*m_r.^2;
        dd_c = Nx*[2*m_d_c(2), m_d_c(2:end)-m_d_c(1:end-1), 2*(Nx*(0-c(end-1).*r(end-1)^2) - m_d_c(end))];

        n_cHat = cHat + dt*RS*Vw*(K*c-cHat);
        c = c + dt*(D*dd_c./r.^2 - RS*Psi*(K*c-cHat));
        cHat = n_cHat;

        c(1) = c(2);
        c(end) = 0;

        t = t + dt;
    end
end