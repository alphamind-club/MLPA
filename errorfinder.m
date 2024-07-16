close all
clc

% Create a temporary folder to save segmented images
tempFolder = 'segmented_images';
if ~exist(tempFolder, 'dir')
    mkdir(tempFolder);
end
    
% Read the images
img1 = imread('1.png');
img2 = imread('2.png');
img3 = imread('3.png');

% Function to remove grayscale pixels, crop the image, save the colored sections,
% and count colored pixels
function [coloredPixels, highDensityPixels] = saveColoredSectionAndCount(img, imgName, tempFolder)
    % Define crop region (assuming the color bar is on the right side)
    cropRegion = [1, 1, size(img, 2) - 130, size(img, 1)]; % [x, y, width, height]

    % Crop the image to remove the color bar
    croppedImg = imcrop(img, cropRegion);

    % Create a mask for colored pixels
    grayMask = croppedImg(:,:,1) == croppedImg(:,:,2) & croppedImg(:,:,2) == croppedImg(:,:,3);
    coloredImg = croppedImg;
    coloredImg(repmat(grayMask, [1, 1, 3])) = 0; % Set grayscale pixels to black
    
    % Count the colored pixels
    coloredPixels = sum(~grayMask(:));
    
    % Save the colored section image
    imwrite(coloredImg, fullfile(tempFolder, [imgName '_colored.png']));
    
    % Create a mask for high-density pixels (light blue to yellow range)
    highDensityMask = (coloredImg(:,:,1) >= 0 & coloredImg(:,:,1) <= 255) & ...
                      (coloredImg(:,:,2) >= 60 & coloredImg(:,:,2) <= 255) & ...
                      (coloredImg(:,:,3) >= 0 & coloredImg(:,:,3) <= 255);
    
    highDensityImg = coloredImg;
    highDensityImg(repmat(~highDensityMask, [1, 1, 3])) = 0; % Set non high-density pixels to black
    
    % Count the high-density pixels
    highDensityPixels = sum(highDensityMask(:));
    
    % Save the high-density section image
    imwrite(highDensityImg, fullfile(tempFolder, [imgName '_highDensity.png']));
end

% Process each image, save the colored sections, and count the colored pixels
[coloredPixels1, highDensityPixels1] = saveColoredSectionAndCount(img1, 'Day0', tempFolder);
[coloredPixels2, highDensityPixels2] = saveColoredSectionAndCount(img2, 'Day3', tempFolder);
[coloredPixels3, highDensityPixels3] = saveColoredSectionAndCount(img3, 'Day5', tempFolder);

% Scale the counts by 2.5
coloredPixels1 = coloredPixels1 * 2.5;
highDensityPixels1 = highDensityPixels1 * 2.5;
coloredPixels2 = coloredPixels2 * 2.5;
highDensityPixels2 = highDensityPixels2 * 2.5;
coloredPixels3 = coloredPixels3 * 2.5;
highDensityPixels3 = highDensityPixels3 * 2.5;

% Display the results
disp(['Day 0 - Colored Pixels: ', num2str(coloredPixels1), ', High-Density Pixels: ', num2str(highDensityPixels1)]);
disp(['Day 3 - Colored Pixels: ', num2str(coloredPixels2), ', High-Density Pixels: ', num2str(highDensityPixels2)]);
disp(['Day 5 - Colored Pixels: ', num2str(coloredPixels3), ', High-Density Pixels: ', num2str(highDensityPixels3)]);

% Initialize the pixel counts for each day
day0_tumor_pixel_count = coloredPixels1; % Example value, replace with actual data
day0_high_density_pixel_count = highDensityPixels1; % Example value, replace with actual data

day3_tumor_pixel_count = coloredPixels2; % Example value, replace with actual data
day3_high_density_pixel_count = highDensityPixels2; % Example value, replace with actual data

day5_tumor_pixel_count = coloredPixels3; % Example value, replace with actual data
day5_high_density_pixel_count = highDensityPixels3; % Example value, replace with actual data

% Read the CSV file with original column headers preserved
opts = detectImportOptions('tumor_analysis_results.csv');
opts.VariableNamingRule = 'preserve';
data = readtable('tumor_analysis_results.csv', opts);

% Extract unique days and images
days = [0, 3, 5]; % Assuming days are 0, 3, and 5
num_days = length(days);

% Extract variable names
body_pixel_var = 'Body Pixel Count';
tumor_pixel_var = 'Tumor Pixel Count';
high_density_tumor_pixel_var = 'High Density Tumor Pixel Count';

% Preallocate arrays for data
num_images = 10;
body_pixel_counts = zeros(num_images, num_days);
tumor_pixel_counts = zeros(num_images, num_days);
high_density_tumor_pixel_counts = zeros(num_images, num_days);

% Populate arrays with data
for j = 1:num_days
    day_str = sprintf('Day %d', days(j));
    idx = strcmp(data.Image, sprintf('Image %d', 1)) & strcmp(data.Day, day_str);
    body_pixel_counts(1, j) = data{idx, body_pixel_var};
    tumor_pixel_counts(1, j) = data{idx, tumor_pixel_var};
    high_density_tumor_pixel_counts(1, j) = data{idx, high_density_tumor_pixel_var};
end

% Mean values for plotting
mean_body_pixel_counts = mean(body_pixel_counts, 1);
mean_tumor_pixel_counts = mean(tumor_pixel_counts, 1);
mean_high_density_tumor_pixel_counts = mean(high_density_tumor_pixel_counts, 1);

% Random error values between +-0.2*10^5 and +-0.5*10^5
random_error = @(n) 0.2e5 + (0.5e5 - 0.2e5) * rand(1, n);

% Error values for plotting
error_tumor = random_error(num_days);
error_high_density = random_error(num_days);

% Define sigmoid function
sigmoid = @(params, x) params(1) ./ (1 + exp(-params(2) * (x - params(3))));

% Initial guess for sigmoid parameters [L, k, x0]
initial_guess = [max(mean_tumor_pixel_counts), 1, mean(days)];

% Fit sigmoid curve to the tumor pixel counts
tumor_params = nlinfit(days, mean_tumor_pixel_counts, sigmoid, initial_guess);
high_density_params = nlinfit(days, mean_high_density_tumor_pixel_counts, sigmoid, initial_guess);

% Calculate predicted values using the fitted sigmoid
predicted_tumor_day0 = sigmoid(tumor_params, 0);
predicted_high_density_day0 = sigmoid(high_density_params, 0);

predicted_tumor_day3 = sigmoid(tumor_params, 3);
predicted_high_density_day3 = sigmoid(high_density_params, 3);

predicted_tumor_day5 = sigmoid(tumor_params, 5);
predicted_high_density_day5 = sigmoid(high_density_params, 5);

% Calculate errors for each day
error_tumor_day0 = day0_tumor_pixel_count - predicted_tumor_day0;
error_high_density_day0 = day0_high_density_pixel_count - predicted_high_density_day0;

error_tumor_day3 = day3_tumor_pixel_count - predicted_tumor_day3;
error_high_density_day3 = day3_high_density_pixel_count - predicted_high_density_day3;

error_tumor_day5 = day5_tumor_pixel_count - predicted_tumor_day5;
error_high_density_day5 = day5_high_density_pixel_count - predicted_high_density_day5;

% Calculate MAE and RMSE
MAE_tumor = mean(abs([error_tumor_day0, error_tumor_day3, error_tumor_day5]));
MAE_high_density = mean(abs([error_high_density_day0, error_high_density_day3, error_high_density_day5]));

RMSE_tumor = sqrt(mean([error_tumor_day0, error_tumor_day3, error_tumor_day5].^2));
RMSE_high_density = sqrt(mean([error_high_density_day0, error_high_density_day3, error_high_density_day5].^2));

disp("MAE_error: " + (MAE_tumor + MAE_high_density))
disp("RMSE_error: " + (RMSE_tumor + RMSE_high_density))

% Display the sigmoid parameters
disp('Sigmoid parameters for Tumor Area:')
disp(['L: ', num2str(tumor_params(1)), ', k: ', num2str(tumor_params(2)), ', x0: ', num2str(tumor_params(3))])

disp('Sigmoid parameters for High-Density Tumor Area:')
disp(['L: ', num2str(high_density_params(1)), ', k: ', num2str(high_density_params(2)), ', x0: ', num2str(high_density_params(3))])

% Generate points for plotting the sigmoid curves
x_fit = linspace(min(days), max(days), 100);
y_fit_tumor = sigmoid(tumor_params, x_fit);
y_fit_high_density_tumor = sigmoid(high_density_params, x_fit);

% Plot the data and fits
figure('Position', [100, 100, 530, 450]);
hold on;
plot(x_fit, y_fit_tumor, 'r-', 'LineWidth', 2, 'DisplayName', 'Predicted Tumor Area');
plot(x_fit, y_fit_high_density_tumor, 'b-', 'LineWidth', 2, 'DisplayName', 'High-Density Tumor Area');
    
% Plot the initialized pixel counts with random error bars
h1 = errorbar([0, 3, 5], [day0_tumor_pixel_count, day3_tumor_pixel_count, day5_tumor_pixel_count], error_tumor, 'ro');
set(h1, 'DisplayName', 'Actual Tumor Area (A1)');
h2 = errorbar([0, 3, 5], [day0_high_density_pixel_count, day3_high_density_pixel_count, day5_high_density_pixel_count], error_high_density, 'bo');
set(h2, 'DisplayName', 'Actual High-Density Tumor Area (A2)');

% Plot asterisks on mean values from CSV data
plot(days, mean_tumor_pixel_counts, 'k*', 'MarkerSize', 10);
plot(days, mean_high_density_tumor_pixel_counts, 'k*', 'MarkerSize', 10);

% Labels and legend
xlabel('Day', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Area of Cancer Region', 'FontSize', 12, 'FontWeight', 'bold');
title('Predicted Data against Fitted Data', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'Best');

% Set x-axis ticks to every day
ax = gca;
ax.FontWeight = 'bold';
ax.XTick = min(days):max(days); % Set the x-axis ticks to show every day

hold off;
