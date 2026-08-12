% Solves the coupled PDE model in a radial geometry

function [C, CH] = radial_coupled_pde_model(D, K, RS, Psi, Vw, Nx, dt, times)
    sz = max(size(times, 1), size(times, 2));

    c = zeros(1, Nx+1);
    cHat = c;

    t = 0;
    reportIndex = 1;

    C = zeros(sz, Nx+1);
    CH = C;

    % Initial condition = 1_{x<1/5}
    c(1:floor(Nx/5)) = 1;

    r = linspace(0, 10, Nx+1);
    m_r = 0.5*(r(2:end)+r(1:end-1));

    while (reportIndex <= sz)
        if(t >= times(reportIndex))
            CH(reportIndex, :) = cHat;
            C(reportIndex, :) = c;
            reportIndex = reportIndex + 1;
        end

        m_d_c = Nx*(c(2:end)-c(1:end-1)).*m_r;
        dd_c = Nx*[2*m_d_c(2), m_d_c(2:end)-m_d_c(1:end-1), -2*m_d_c(end)];

        n_cHat = cHat + dt*RS*Vw*(K*c-cHat);
        c = c + dt*(D*dd_c./r - RS*Psi*(K*c-cHat));
        cHat = n_cHat;

        c(1) = c(2);
        c(end) = c(end-1);

        t = t + dt;
    end
end