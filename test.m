close all

% Create a random growth rate map
rng('shuffle'); % For reproducibility
sigma = 1; % Standard deviation of the Gaussian filter
rawMap = imgaussfilt(rand(gridSize, gridSize), sigma); % Smooth the random values

% Normalize the map to the range [0, 1]
minVal = min(rawMap(:));
maxVal = max(rawMap(:));
polarizedMap = ((rawMap - minVal) / (maxVal - minVal)).^2;

% Transform the values based on the specified thresholds
growthRateMap = zeros(size(polarizedMap)); % Initialize the growth rate map

% Apply the thresholds
growthRateMap(polarizedMap > 0.6) = 1;      % Values above 0.6 set to 1
growthRateMap(polarizedMap > 0.3 & polarizedMap <= 0.6) = 0.5; % Values between 0.3 and 0.6 set to 0.5
growthRateMap(polarizedMap <= 0.3) = 0;     % Values below or equal to 0.3 set to 0

% Display the growth rate map
figure('Name', 'Polarized Map', 'NumberTitle', 'off');
imagesc(polarizedMap);
colorbar;
title('Polarized Map');
xlabel('X');
ylabel('Y');

figure('Name', 'Growth Rate Map', 'NumberTitle', 'off');
imagesc(growthRateMap);
colorbar;
title('Growth Rate Map');
xlabel('X');
ylabel('Y');
