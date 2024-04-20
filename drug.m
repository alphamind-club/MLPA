close all;

% Parameters
rise_rate = 2;                   % Rate of rise
decay_rate = 0.05;               % Decay rate
num_spikes = 3;                  % Number of spikes

% Time vector
t = 0:0.1:150;                   % Time points from 0 to 50 with a step size of 0.1

% Initialize arrays for spikes
spike_amplitude = [50, 50, 50];  % Amplitude of each spike
spike_time = [10, 60, 110];        % Time at which each spike occurs

% Initialize drug concentration
drug_concentration = zeros(size(t));

% Add each spike to the drug concentration
for i = 1:num_spikes
    gradual_rise = 1 ./ (1 + exp(-rise_rate * (t - spike_time(i))));
    drug_concentration = drug_concentration + gradual_rise .* spike_amplitude(i) .* exp(-decay_rate * (t - spike_time(i)));
end

% Plot the results
figure;
plot(t, drug_concentration, 'LineWidth', 2);
title('Drug Metabolization Graph');
xlabel('Time');
ylabel('Drug Concentration');
grid on;