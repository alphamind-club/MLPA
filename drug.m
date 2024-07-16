function tumor_growth_gui()
    clc
    close all

    % Create figure
    fig = figure('Position', [200, 200, 400, 300], 'Name', 'Tumor Growth Model GUI');
    
    % Initialize parameters
    rise_rate = 10;
    decay_rate = 0.5;
    initialSize = 60;  % Adjusted initial size for the new scale
    growthRate = 0.025;  % Adjusted growth rate for the new scale
    carryingCapacity = 12000;
    molecularWeight = 500;
    spike_amplitude = 35;
    spike_time = [5, 9, 13, 17];
    
    % Create UI controls
    uicontrol('Style', 'text', 'String', 'Rise Rate:', 'Position', [20, 250, 100, 20]);
    riseRateEdit = uicontrol('Style', 'edit', 'String', num2str(rise_rate), 'Position', [150, 250, 100, 20]);
    
    uicontrol('Style', 'text', 'String', 'Decay Rate:', 'Position', [20, 220, 100, 20]);
    decayRateEdit = uicontrol('Style', 'edit', 'String', num2str(decay_rate), 'Position', [150, 220, 100, 20]);
    
    uicontrol('Style', 'text', 'String', 'Initial Tumor Size:', 'Position', [20, 190, 120, 20]);
    initialSizeEdit = uicontrol('Style', 'edit', 'String', num2str(initialSize), 'Position', [150, 190, 100, 20]);
    
    uicontrol('Style', 'text', 'String', 'Tumor Growth Rate:', 'Position', [20, 160, 130, 20]);
    growthRateEdit = uicontrol('Style', 'edit', 'String', num2str(growthRate), 'Position', [150, 160, 100, 20]);
    
    uicontrol('Style', 'text', 'String', 'Carrying Capacity:', 'Position', [20, 130, 120, 20]);
    carryingCapacityEdit = uicontrol('Style', 'edit', 'String', num2str(carryingCapacity), 'Position', [150, 130, 100, 20]);
    
    uicontrol('Style', 'text', 'String', 'Molecular Weight:', 'Position', [20, 100, 120, 20]);
    molecularWeightEdit = uicontrol('Style', 'edit', 'String', num2str(molecularWeight), 'Position', [150, 100, 100, 20]);
    
    uicontrol('Style', 'text', 'String', 'Spike Amplitude:', 'Position', [20, 70, 120, 20]);
    spikeAmplitudeEdit = uicontrol('Style', 'edit', 'String', num2str(spike_amplitude), 'Position', [150, 70, 100, 20]);
    
    uicontrol('Style', 'text', 'String', 'Spike Time:', 'Position', [20, 40, 120, 20]);
    spikeTimeEdit = uicontrol('Style', 'edit', 'String', num2str(spike_time), 'Position', [150, 40, 100, 20]);
    
    uicontrol('Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [150, 10, 100, 20], 'Callback', @runSimulation);
    
    function runSimulation(~, ~)
        % Get parameter values
        rise_rate = str2double(get(riseRateEdit, 'String'));
        decay_rate = str2double(get(decayRateEdit, 'String'));
        initialSize = str2double(get(initialSizeEdit, 'String'));
        growthRate = str2double(get(growthRateEdit, 'String'));
        carryingCapacity = str2double(get(carryingCapacityEdit, 'String'));
        molecularWeight = str2double(get(molecularWeightEdit, 'String'));
        spike_amplitude = str2double(get(spikeAmplitudeEdit, 'String'));
        spike_time = str2num(get(spikeTimeEdit, 'String'));  %#ok<ST2NM>
        
        % Check for empty inputs and retain old parameters
        if isempty(rise_rate)
            rise_rate = 10;
        end
        if isempty(decay_rate)
            decay_rate = 1;
        end
        if isempty(initialSize)
            initialSize = 60;
        end
        if isempty(growthRate)
            growthRate = 0.03;
        end
        if isempty(carryingCapacity)
            carryingCapacity = 6000;
        end
        if isempty(molecularWeight)
            molecularWeight = 500;
        end
        if isempty(spike_amplitude)
            spike_amplitude = 35;
        end
        if isempty(spike_time)
            spike_time = [2, 7, 12, 17];
        end
        
        % Run simulation with the updated parameters
        simulateTumorGrowth(rise_rate, decay_rate, initialSize, growthRate, carryingCapacity, molecularWeight, spike_amplitude, spike_time);
    end
end

function simulateTumorGrowth(rise_rate, decay_rate, initialSize, growthRate, carryingCapacity, molecularWeight, spike_amplitude, spike_time)
    % Simulation parameters
    timesteps = 300;  % Number of iterations
    step_size = 0.1;
    t_real = 0:step_size:20;  % Time axis from 0 to 20 (real time units)
    t_iter = linspace(0, timesteps, length(t_real));  % Corresponding iteration axis from 0 to 300
    
    num_spikes = length(spike_time);

    % Initialize drug concentration
    drug_concentration = zeros(size(t_real));
    noise_amplitude = 0;  % Noise amplitude for drug concentration
    
    % Add each spike to the drug concentration with noise
    for i = 1:num_spikes
        gradual_rise = 1 ./ (1 + exp(-rise_rate * (t_real - spike_time(i))));
        drug_concentration = drug_concentration + gradual_rise .* spike_amplitude .* exp(-decay_rate * (t_real - spike_time(i))) / molecularWeight;
    end
    
    % Tumor growth model (Logistic growth equation)
    tumorSize = zeros(1, length(t_real));
    tumorSize(1) = initialSize;

    % Tumor growth without treatment
    tumorSizeNoTreatment = zeros(1, length(t_real));
    tumorSizeNoTreatment(1) = initialSize;

    for i = 2:length(t_real)
        if i >= 75
            noise_amplitude = 5;
        end
        noise = noise_amplitude * randn;  % Additive Gaussian noise for each time step
        % Update tumor size with noise
        drug_effect = drug_concentration(i-1) * 5700;
        tumorSize(i) = tumorSize(i-1) + growthRate * (tumorSize(i-1) - drug_effect) * (1 - tumorSize(i-1) / carryingCapacity) + noise;
        
        % Update tumor size without treatment
        tumorSizeNoTreatment(i) = tumorSizeNoTreatment(i-1) + growthRate * tumorSizeNoTreatment(i-1) * (1 - tumorSizeNoTreatment(i-1) / carryingCapacity) + noise;

        if tumorSize(i) <= 0
            tumorSize(i) = 0;
        end
        if tumorSizeNoTreatment(i) <= 0
            tumorSizeNoTreatment(i) = 0;
        end
    end
    
    % Plot drug concentration
    figure('Position', [200, 200, 600, 200]);  % Half the height
    plot(t_iter, drug_concentration, 'r-', 'LineWidth', 2);
    ylabel('Drug Concentration', 'FontSize', 12, 'FontWeight', 'bold');
    title('Drug Metabolization', 'FontSize', 14, 'FontWeight', 'bold');
    ylim([0 0.07]);

    % Plot tumor growth with and without treatment
    figure('Position', [200, 200, 600, 200]);  % Half the height
    plot(t_iter, tumorSizeNoTreatment, 'b-', 'LineWidth', 2); % Dotted line for no treatment
    hold on;
    plot(t_iter, tumorSize, 'r-', 'LineWidth', 2);
    hold off;
    legend('With Treatment', 'Without Treatment');
    ylabel('Tumor Size', 'FontSize', 12, 'FontWeight', 'bold');
    title('Tumor Growth with Treatment', 'FontSize', 14, 'FontWeight', 'bold');
    ylim([0 6000]);  % Set y-axis limits from 0 to 6000 tumor cells
end
