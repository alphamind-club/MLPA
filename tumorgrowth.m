% Tumor Growth Model

% Parameters
initialSize = 1;    % Initial tumor size
growthRate = 0.1;   % Tumor growth rate
carryingCapacity = 100;  % Maximum tumor size

% Time vector
timesteps = 150;
t = 0:1:(timesteps-1);  % Adjusted to have the same length as timesteps

% Tumor growth model (Logistic growth equation)
tumorSize = zeros(1, timesteps);
tumorSize(1) = initialSize;

for i = 2:timesteps
    tumorSize(i) = tumorSize(i-1) + growthRate * tumorSize(i-1) * (1 - tumorSize(i-1)/carryingCapacity);
end

% Plotting the results
figure;
plot(t, tumorSize, 'LineWidth', 2);
title('Tumor Growth Model');
xlabel('Time');
ylabel('Tumor Size');
grid on;