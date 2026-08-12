K = 10;
N = 3;

func_form = @(t, deff, f)(2.*deff.*(sin(pi*f/2) - sin(3*pi*f/2).*exp(-2*deff*pi^2*t)).*exp(-deff*pi^2*t/4));

total_flux = @(t, C)(sum(C(1:N).*func_form(t, C(N+1:2*N), C(2*N+1:end))));

times = linspace(1e-3, 20, 500);

coeffs = [1, 10, 100, 1, 1/10, 1/100, 0.55, 0.5, 0.45];
plot(times, arrayfun(@(t)total_flux(t, coeffs), times));
hold on;
plot(times, coeffs(1)*func_form(times, coeffs(4), coeffs(7)));
plot(times, coeffs(2)*func_form(times, coeffs(5), coeffs(8)));
plot(times, coeffs(3)*func_form(times, coeffs(6), coeffs(9)));

clf;

coeffs = [12, 14, 77, 0.05, 0.5, 0.02, 0.3, 0.2, 0.3];
plot(times, arrayfun(@(t)total_flux(t, coeffs), times));
hold on;
plot(times, coeffs(1)*func_form(times, coeffs(4), coeffs(7)));
plot(times, coeffs(2)*func_form(times, coeffs(5), coeffs(8)));
plot(times, coeffs(3)*func_form(times, coeffs(6), coeffs(9)));