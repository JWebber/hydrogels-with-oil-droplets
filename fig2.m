starts = [150, 150] + [0, 0; -1, 0; -1, -1; 0, -1; 1, 0; 0, 1; 1, 1; 1, -1; -1, 1];
starts = repmat(starts, [100, 1]);
N = 2500;
repunits = 30;

tiledlayout(3, 5, 'TileSpacing','compact','Padding','tight');

i = 1;
for K = [0.1 1 10 100 1000]
    nexttile(i)
    paths = agent_based_diffusion(starts, N, K, 10, 2, repunits);
    rectangle('Position', [0, 0, 10*repunits, 10*repunits], 'edgecolor', 'black', 'linewidth', 0.5);
    hold on;

    for j = 1 : size(starts, 1)
        plot(paths(:, j, 1), paths(:, j, 2));
        hold on;
    end
    axis off;
    xlim([0 10*repunits]);
    ylim([0 10*repunits]);
    drawnow;

    exportgraphics(gca, strcat('agent_plots/paths_', num2str(K), '.png'), 'Resolution', 800);
    cla;

    pause(0.01);

    radii = sqrt((squeeze(paths(:, :, 1))-paths(1, :, 1)).^2 + (squeeze(paths(:, :, 2))-paths(1, :, 2)).^2);
    mr = mean(radii, 2, "omitnan");
    stdr = std(radii, 0, 2, "omitnan");
    
    nexttile(i+5);

    [x, y] = minplot(1:N, mean(radii, 2, "omitnan"), 1e-3);
    loglog(x, y);
    matlab2tikz(strcat('agent_plots/mean_radii_', num2str(K), '.tikzdata'));

    cla;

    i = i+1;
end