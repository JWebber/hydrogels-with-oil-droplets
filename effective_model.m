% Models 1d diffusion for figure 5

function [C, F] = effective_model(D, fillfrac, K, RS, Psi, Vw, Nx, dt, times)
    sz = max(size(times, 1), size(times, 2));

    c = [ones(1, floor(Nx*fillfrac)) zeros(1, Nx+1-floor(Nx*fillfrac))];

    t = 0;
    reportIndex = 1;

    C = zeros(sz, Nx+1);
    F = zeros(sz, 1);

    Deff = (1-K*Psi/(Vw+K*Psi))*D;

    while (reportIndex <= sz)
        if(t >= times(reportIndex))
            C(reportIndex, :) = c;
            F(reportIndex) = -Nx*Deff*(c(end)-c(end-1));
            reportIndex = reportIndex + 1;
        end

        m_d_c = Nx*(c(2:end)-c(1:end-1));
        dd_c = Nx*[2*m_d_c(2), m_d_c(2:end)-m_d_c(1:end-1), 2*(Nx*(0-c(end-1))-m_d_c(end))];

        c = c + dt*Deff*dd_c;

        c(1) = c(2);
        c(end) = 0;

        t = t + dt;
    end
end