starts = [150, 150] + [0, 0; -1, 0; -1, -1; 0, -1; 1, 0; 0, 1; 1, 1; 1, -1; -1, 1];
starts = repmat(starts, [50, 1]);
N = 2500;
repunits = 30;

Ks = linspace(0.01, 100, 100);

radii = zeros(100, N);
fits = zeros(100, 1);

for i = 1 : 100
    paths = agent_based_diffusion(starts, N, Ks(i), 10, 2, repunits);

    c_radii = sqrt((squeeze(paths(:, :, 1))-paths(1, :, 1)).^2 + (squeeze(paths(:, :, 2))-paths(1, :, 2)).^2);
    mr = mean(c_radii, 2, "omitnan");

    prevG = 1;
    newG = 2;

    stp = 1;

    while(abs(newG-prevG)/prevG > 0.01 && stp < N - 200)
        fitfun = @(r)(sqrt(r)*sqrt(stp:N) - mr(stp:end));

        prevG = newG;
        newG = lsqnonlin(fitfun, prevG);

        stp = stp+100;
    end

    radii(i, :) = mr;
    fits(i) = newG;

    plot(Ks(i), newG, 'r.');
    hold on;
    pause(0.01);

    save('fit_diffs.mat', "radii", "fits");
end