function tumor_growth_gui()
    clc
    close all

    % Create figure
    fig = figure('Position', [200, 200, 400, 300], 'Name', 'Tumor Growth Model GUI');
    
    % Initialize parameters
    rise_rate = 10;
    decay_rate = 0.5;
    initialSize = 0.01;
    growthRate = 0.0235;
    carryingCapacity = 250;
    molecularWeight = 500;
    spike_amplitude = 35;
    spike_time = [5, 10, 15];
    
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
            initialSize = 0.01;
        end
        if isempty(growthRate)
            growthRate = 0.0235;
        end
        if isempty(carryingCapacity)
            carryingCapacity = 250;
        end
        if isempty(molecularWeight)
            molecularWeight = 500;
        end
        if isempty(spike_amplitude)
            spike_amplitude = 35;
        end
        if isempty(spike_time)
            spike_time = [5, 10, 15];
        end
        
        % Run simulation with the updated parameters
        simulateTumorGrowth(rise_rate, decay_rate, initialSize, growthRate, carryingCapacity, molecularWeight, spike_amplitude, spike_time);
    end
end

function simulateTumorGrowth(rise_rate, decay_rate, initialSize, growthRate, carryingCapacity, molecularWeight, spike_amplitude, spike_time)
    % Simulation parameters
    timesteps = 200;
    step_size = 0.1;
    t = 0:step_size:(timesteps-1)*step_size;
    num_spikes = length(spike_time);
    
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
end