i = 1;
bg = {"#cccce6", "#9999cc", "#6666b2", "#333399", "#000080"};

mean_radii = NaN*zeros(2500, 6);
mean_radii(:, 1) = [1:2500]';

for K = [0 1 10 100 1000]
    clf;
    [~, allpos, means] = agent_sim(K, 500, 2, 10, 21, 2500, 50);
    rectangle('position', [-120, -120, 240, 240], 'edgecolor', 'none', 'facecolor', bg{i});
    hold on;
    plot(allpos(:, :, 1), allpos(:, :, 2));
    xlim([-120 120]);
    ylim([-120 120]);
    drawnow;
    axis off;
    exportgraphics(gca, strcat('K_', num2str(K), '.png'), 'Resolution', 500);
    i = i+1;

    mean_radii(:, i) = (means - means(1))';
end

writematrix(mean_radii, "mean_radii.csv");