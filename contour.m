close all

% Data preparation
tumorChance = [0.055, 0.04, 0.06, 0.045, 0.055, 0.05, 0.04, 0.06, 0.045, 0.05, 0.055, 0.055, 0.04, 0.04, 0.045, 0.05, 0.06, 0.06, 0.045, 0.05];
metastasisChance = [0.0002, 0.0004, 0.0001, 0.0003, 0.0004, 0.0002, 0.0001, 0.0004, 0.0002, 0.0003, 0.0003, 0.0001, 0.0002, 0.0003, 0.0001, 0.0001, 0.0002, 0.0003, 0.0004, 0.0004];
MAE_Error = [77469.3667, 257688.5333, 109317.0895, 251457.7, 142937.7105, 142151.8667, 304775.2, 161231.132, 278241.731, 122359.442, 108214.51, 142315.21, 301451.214, 271953.6345, 265513.723, 218245.796, 127623.416, 166514.232, 221045.929, 162953.14];

%tumorChance = [0.045, 0.055, 0.035, 0.03, 0.065, 0.03, 0.065];
%metastasisChance = [0.0003, 0.0002, 0.0004, 0.0002, 0.0004, 0.0004, 0.0002];
%MAE_Error = [48318, 141562, 192638, 275162, 282715, 277390, 251842];

% Create a grid for contour plot
[tumorGrid, metastasisGrid] = meshgrid(linspace(min(tumorChance), max(tumorChance), 100), linspace(min(metastasisChance), max(metastasisChance), 100));

% Interpolate the MAE_Error values over the grid
MAE_Error_Grid = griddata(tumorChance, metastasisChance, MAE_Error, tumorGrid, metastasisGrid, 'cubic');

% Plot the contour plot
figure('Position', [100, 100, 400, 300]); % Adjust the figure size
contourf(tumorGrid, metastasisGrid, MAE_Error_Grid, 20);
colorbar;
hold on;

% Plot the data points
% scatter(tumorChance, metastasisChance, 100, MAE_Error, 'filled', 'MarkerEdgeColor', 'k');
xlabel('Pgrowth', 'FontName', 'Helvetica', 'FontWeight', 'bold');
ylabel('Pmetastasis', 'FontName', 'Helvetica', 'FontWeight', 'bold');
title('Contour Plot of Model Error', 'FontName', 'Helvetica');
set(gca, 'FontName', 'Helvetica'); % Set the font for the axes
hold off;
