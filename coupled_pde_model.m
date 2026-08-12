% Models 1D diffusion for figure 5

function [C, F] = coupled_pde_model(D, K, RS, Psi, Vw, Nx, dt, times)
    sz = max(size(times, 1), size(times, 2));

    c = [ones(1, floor(Nx/2)) zeros(1, Nx+1-floor(Nx/2))];
    cHat = c;

    t = 0;
    reportIndex = 1;

    C = zeros(sz, Nx+1);
    F = zeros(sz, 1);

    while (reportIndex <= sz)
        if(t >= times(reportIndex))
            C(reportIndex, :) = c;
            F(reportIndex) = -Nx*D*(c(end)-c(end-1));
            reportIndex = reportIndex + 1;
        end

        m_d_c = Nx*(c(2:end)-c(1:end-1));
        dd_c = Nx*[2*m_d_c(2), m_d_c(2:end)-m_d_c(1:end-1), 2*(Nx*(0-c(end-1))-m_d_c(end))];

        n_cHat = cHat + dt*RS*Vw*(K*c-cHat);
        c = c + dt*(D*dd_c - RS*Psi*(K*c-cHat));
        cHat = n_cHat;

        c(1) = c(2);
        c(end) = 0;

        t = t + dt;
    end
end