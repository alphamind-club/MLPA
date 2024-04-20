clc
close all

% Parameters
rise_rate = 10;                    % Rise rate of drug
decay_rate = 0.5;                  % Decay rate of drug

initialSize = 0.01;                % Initial tumor size in cm
growthRate = 0.245;                % Tumor growth rate
carryingCapacity = 250;            % Maximum tumor size
molecularWeight = 500;             % Molecular weight of the drug in g/mol

% Time vector
timesteps = 20;
t = 0:1:(timesteps-1);                   % Time points from 0 to 50 with a step size of 0.1

% Initialize arrays for spikes
spike_amplitude = 35;  % Amplitude of each spike
spike_time = [5,10,15];      % Time at which each spike occurs
num_spikes = length(spike_time);                  % Number of spikes

% Initialize drug concentration
drug_concentration = zeros(size(t));

% Add each spike to the drug concentration
for i = 1:num_spikes
    gradual_rise = 1 ./ (1 + exp(-rise_rate * (t - spike_time(i))));
    drug_concentration = drug_concentration + gradual_rise .* spike_amplitude .* exp(-decay_rate * (t - spike_time(i))) / molecularWeight;
end

% Tumor growth model (Logistic growth equation)
tumorSize = zeros(1, timesteps);
tumorSize(1) = initialSize;

for i = 2:timesteps
    tumorSize(i) = tumorSize(i-1) + growthRate * (tumorSize(i-1) - drug_concentration(i-1)) * (1 - tumorSize(i-1) / carryingCapacity);
    
    if tumorSize(i) <= 0
        tumorSize(i) = 0;
    end
end

% Plot drug concentration
figure;
plot(t, drug_concentration, 'LineWidth', 2);
title('Drug Metabolization Graph');
xlabel('Time');
ylabel('Drug Concentration in the Blood');
grid on;

% Plot tumor growth
figure;
plot(t, tumorSize, 'LineWidth', 2);
title('Tumor Growth Model');
xlabel('Time');
ylabel('Tumor Size');
grid on;

% Suggested dosage rate calculation
amplitude_test = spike_amplitude;
while tumorSize(timesteps) > 0
    drug_concentration = zeros(size(t));
    tumorSize = zeros(1, timesteps);
    tumorSize(1) = initialSize;

    amplitude_test = amplitude_test + 1;

    for i = 1:num_spikes
        gradual_rise = 1 ./ (1 + exp(-rise_rate * (t - spike_time(i))));
        drug_concentration = drug_concentration + gradual_rise .* amplitude_test .* exp(-decay_rate * (t - spike_time(i))) / molecularWeight;
    end

    for i = 2:timesteps
        tumorSize(i) = tumorSize(i-1) + growthRate * (tumorSize(i-1) - drug_concentration(i-1)) * (1 - tumorSize(i-1) / carryingCapacity);
    
        if tumorSize(i) <= 0
            tumorSize(i) = 0;
        end
    end
end

disp("The suggested dosage rate should be " + amplitude_test + " to prevent relapse.")
