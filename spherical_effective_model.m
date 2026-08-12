function [C, F] = spherical_effective_model(D, K, RS, Psi, Vw, Nx, dt, times)
    sz = max(size(times, 1), size(times, 2));

    c = ones(1, Nx+1);

    t = 0;
    reportIndex = 1;

    C = zeros(sz, Nx+1);
    F = zeros(sz, 1);

    Deff = (1-K*Psi/(Vw+K*Psi))*D;

    r = linspace(0, 10, Nx+1);
    m_r = 0.5*(r(2:end)+r(1:end-1));

    while (reportIndex <= sz)
        if(t >= times(reportIndex))
            C(reportIndex, :) = c;
            F(reportIndex) = -Nx*Deff*(c(end)-c(end-1));
            reportIndex = reportIndex + 1;
        end

        m_d_c = Nx*(c(2:end)-c(1:end-1)).*m_r.^2;
        dd_c = Nx*[2*m_d_c(2), m_d_c(2:end)-m_d_c(1:end-1), 2*(Nx*(0-c(end-1).*r(end-1)^2) - m_d_c(end))];

        c = c + dt*Deff*dd_c./r.^2;

        c(1) = c(2);
        c(end) = 0;

        t = t + dt;
    end
end