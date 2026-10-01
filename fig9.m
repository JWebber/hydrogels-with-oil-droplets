R = 2.821; % 25% volume effect
Deffs = zeros(10,2);

for j = 1 : 10
    [Deffs(j, :), ~, ~] = agent_sim(1, 2500, R, 10, 45, 3000, 2500);
end

De1 = sum(Deffs(:, 1).*Deffs(:, 2))/sum(Deffs(:, 2));

diffusivities = zeros(50, 3);
diffusivities(:, 1) = linspace(0, 1, 50)';

ratio = R^2*pi/10^2/(1-R^2*pi/10^2);

for i = 1 : 50
    Deffs = zeros(10,2);

    for j = 1 : 10
        [Deffs(j, :), ~, ~] = agent_sim(diffusivities(i, 1), 2500, R, 10, 45, 3000, 2500);
    end

    effective_diff = sum(Deffs(:, 1).*Deffs(:, 2))/sum(Deffs(:, 2));

    plot(diffusivities(i, 1), effective_diff/De1, 'b.');
    diffusivities(i, 2) = effective_diff/De1;

    hold on;
    plot(diffusivities(i, 1), (1/(1+ratio*diffusivities(i, 1)))/(1/(1+ratio)), 'r.');
    diffusivities(i, 3) = (1/(1+ratio*diffusivities(i, 1)))/(1/(1+ratio));
    set(gca, 'xscale', 'log');
    drawnow;
end

writematrix(diffusivities, strcat('R_', num2str(R), '.csv'));