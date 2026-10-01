Rs = [0.3989, 0.5, 0.5642, 1, 1.2616, 1.7841];

Z = zeros(6, 100);
X = zeros(6, 100);
Y = zeros(6, 100);

for i = 1 : 6
    matval = readmatrix(strcat("R_", num2str(Rs(i)), ".csv"));

    errs = abs(matval(:, 2)- matval(:, 3))./matval(:, 3);

    Z(i, :) = interp1(matval(:, 1)', errs', linspace(1, 75, 100));
    X(i, :) = pi*Rs(i)^2/100;
    Y(i, :) = linspace(1, 75, 100);
end

contour(X, Y, Z, [0.05, 0.1, 0.25, 0.5, 1], 'LineWidth', 2, 'EdgeColor','red');
matlab2tikz('error_contours.tikzdata');