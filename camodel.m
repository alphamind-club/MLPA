function camodel(breakThreshold, tumorChance, deathChance, mutationChance, metastasisChance)
    % Close all existing figures
    close all;

    % Constant parameters
    gridSize = 200;
    initialTumorSize = 2;
    neighborhoodSize = 1;
    growhoodSize = 3;
    vesselHumanizer = 0.5;
    growChance = 1;
    minTumorCellCountForMetastasis = 50;
    treatment_time = [];

    % Initialize grid
    grid = ones(gridSize, gridSize);

    % Initialize tumor
    center = [floor(gridSize / 2), floor(gridSize / 2)];
    [x, y] = meshgrid(1:gridSize, 1:gridSize);
    circleMask = (x - center(1)).^2 + (y - center(2)).^2 <= (initialTumorSize / 2)^2;
    grid(circleMask) = 2;

    % Create blood vessel - rectangular shape starting from the edge
    vesselWidth = 5;
    vesselHeight = 150;
    startRow = 1;
    startCol = 140;
    vesselMask = false(gridSize, gridSize);
    vesselMask(startRow:(startRow + vesselHeight - 1), startCol:(startCol + vesselWidth - 1)) = true;
    grid(vesselMask) = 5; % Assign a new value (e.g., 5) for blood vessel

    % Create a random growth rate map
    rng(1234); % For reproducibility
    sigma = 1; % Standard deviation of the Gaussian filter
    rawMap = imgaussfilt(rand(gridSize, gridSize), sigma); % Smooth the random values

    % Normalize the map to the range [0, 1]
    minVal = min(rawMap(:));
    maxVal = max(rawMap(:));
    polarizedMap = ((rawMap - minVal) / (maxVal - minVal)).^2;

    % Transform the values based on the specified thresholds
    growthRateMap = zeros(size(polarizedMap)); % Initialize the growth rate map

    % Apply the thresholds
    growthRateMap(polarizedMap > 0.6) = 0.02;      % Values above 0.6 set to 1
    growthRateMap(polarizedMap > 0.3 & polarizedMap <= 0.6) = 0.02; % Values between 0.3 and 0.6 set to 0.5
    growthRateMap(polarizedMap <= 0.3) = 1;     % Values below or equal to 0.3 set to 0

    % Adjust tumor chance based on the growth rate map
    adjustedTumorChance = tumorChance * growthRateMap;

    % GUI setup
    hFig = figure('Name', 'Cancer Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);
    hAx = axes('Parent', hFig, 'Position', [0.05, 0.1, 0.4, 0.8]);
    hText = uicontrol('Style', 'text', 'Parent', hFig, 'Position', [700, 300, 600, 200], 'FontSize', 12, 'HorizontalAlignment', 'left');

    % Display initial parameters
    paramText = sprintf(['Simulation Parameters:\n', ...
        'Break Threshold: %f\n', ...
        'Tumor Chance: %f\n', ...
        'Death Chance: %f\n', ...
        'Mutation Chance: %f\n', ...
        'Metastasis Chance: %f\n'], breakThreshold, tumorChance, deathChance, mutationChance, metastasisChance);
    set(hText, 'String', paramText);

    % Visualization setup
    colormap([1 1 1; 0 0 1; 1 0.75 0.7; 0.5 0 0; 1 0.75 0.7; 0.5 0.5 1]); % 2 is tumor, 3 is grown vessel, 5 is original vessel

    % Simulation parameters
    numIterations = 300;
    numberOfBreaks = 0;
    normalTumorCount = zeros(1, numIterations);
    mutatedTumorCount = zeros(1, numIterations);
    netSurvival = zeros(1, numIterations); % Initialize net survival array
    maxDistance = abs(center(2) - startCol);
    metastasis = 0;

    % Initialize density map
    densityMap = zeros(gridSize, gridSize);

    % Define neighborhood offsets
    [dx, dy] = meshgrid(-neighborhoodSize:neighborhoodSize, -neighborhoodSize:neighborhoodSize);
    neighborhoodOffsets = [dx(:), dy(:)];
    neighborhoodOffsets(all(neighborhoodOffsets == 0, 2), :) = [];

    % Define growhood offsets
    [gx, gy] = meshgrid(-growhoodSize:growhoodSize, -growhoodSize:growhoodSize);
    growhoodOffsets = [gx(:), gy(:)];
    growhoodOffsets(all(growhoodOffsets == 0, 2), :) = [];

    % Run simulation
    for iteration = 1:numIterations
        % Visualization
        axes(hAx);
        image(grid);
        set(gca, 'Color', [1 1 1]); % Set background to white
        title(['Iteration: ' num2str(iteration)]);
        drawnow;

        % Update grid
        newGrid = grid;

        % Tumor growth
        growMask = false(size(grid));
        for i = 1:gridSize
            for j = 1:gridSize
                if grid(i, j) == 2 || grid(i, j) == 6
                    % Update density map
                    densityMap(i, j) = densityMap(i, j) + 1;

                    neighbors = bsxfun(@plus, neighborhoodOffsets, [i, j]);
                    validNeighbors = neighbors(all(neighbors > 0 & neighbors <= gridSize, 2), :);
                    for n = 1:size(validNeighbors, 1)
                        ni = validNeighbors(n, 1);
                        nj = validNeighbors(n, 2);
                        if grid(ni, nj) == 1
                            % Check if adjacent to blood vessel
                            adjacentToBloodVessel = any(grid(sub2ind(size(grid), validNeighbors(:,1), validNeighbors(:,2))) == 3 | grid(sub2ind(size(grid), validNeighbors(:,1), validNeighbors(:,2))) == 5);
                            if adjacentToBloodVessel
                                growthChance = adjustedTumorChance(ni, nj) * 1.5; % Increase growth rate near blood vessels
                            else
                                growthChance = adjustedTumorChance(ni, nj);
                            end
                            if rand <= growthChance
                                growMask(ni, nj) = true;
                            end
                        end
                    end
                end
                
                % Create a break
                neighborhood = grid(max(1, i-neighborhoodSize):min(gridSize, i+neighborhoodSize), max(1, j-neighborhoodSize):min(gridSize, j+neighborhoodSize));
                emptySpace = sum(neighborhood(:) == 1);
                tumorCount = sum(neighborhood(:) == 2) + sum(neighborhood(:) == 6);
                if grid(i,j) == 5 && emptySpace ~= 0 && (i < 160 && i > 40)
                    if numberOfBreaks == 0 && rand <= breakThreshold
                        newGrid(i,j) = 4;
                        numberOfBreaks = numberOfBreaks + 1;
                    elseif rand <= breakThreshold / (10 * numberOfBreaks)
                        newGrid(i,j) = 4;
                        numberOfBreaks = numberOfBreaks + 1;
                    end
                end

                % Select the closest tumor cell to grow to
                tumorDistance = 1e+30;
                if grid(i,j) == 4 && tumorCount == 0
                    for a = 1:gridSize
                        for b = 1:gridSize
                            if grid(a,b) == 2 || grid(a,b) == 6
                                tempDistance = sqrt((a - i)^2 + (b - j)^2);
                                if tempDistance < tumorDistance
                                    tumorA = a;
                                    tumorB = b;
                                    tumorDistance = tempDistance;
                                end
                            end
                        end
                    end
                end

                % Grow towards the closest tumor cell
                growDistance = 1e+30;
                growC = 0;
                growD = 0;
                if grid(i,j) == 4 && rand <= growChance && tumorCount == 0
                    for c = i-1:i+1
                        for d = j-1:j+1
                            tempGrowDistance = sqrt((c - tumorA)^2 + (d - tumorB)^2);
                            if tempGrowDistance < growDistance
                                growDistance = tempGrowDistance;
                                growC = c;
                                growD = d;
                            end
                        end
                    end
                    if rand <= vesselHumanizer
                        growC = i + randi([-1, 1]);
                        growD = j + randi([-1, 1]);
                    end

                    bloodVesselWidth = floor(5 * abs(j - center(2)) / maxDistance);
                    for w = -bloodVesselWidth:bloodVesselWidth
                        newGrid(i + w, j) = 3;
                    end

                    % Make sure it doesn't grow outside the box
                    if growC < 1
                        growC = 1;
                    elseif growC > 200
                        growC = 200;
                    end
                    if growD < 1
                        growD = 1;
                    elseif growD > 200
                        growD = 200;
                    end

                    newGrid(growC, growD) = 4;
                end

                % Stop if touch tumor
                if grid(i,j) == 4 && tumorCount ~= 0
                    newGrid(i, j) = 3;
                    metastasis = 1;
                end
            end
        end
        newGrid(growMask) = 2;

        % Mutate
        mutateMask = grid == 2 & rand(gridSize, gridSize) <= mutationChance;
        newGrid(mutateMask) = 6;

        % Drug administration
        if ismember(iteration, treatment_time)
            deathMask = (grid == 2 & rand(gridSize, gridSize) <= deathChance) | ...
                        (grid == 6 & rand(gridSize, gridSize) <= deathChance / 2);
            newGrid(deathMask) = 1;
        end

        % Metastasis
        if sum(grid(:) == 2) + sum(grid(:) == 6) >= minTumorCellCountForMetastasis && metastasis == 1
            numMetastasisSites = randi(3); % Random number of metastasis sites (between 1 and 3)
            for n = 1:numMetastasisSites
                [metaRow, metaCol] = find(grid == 5 | grid == 3); % Find blood vessel cells
                adjacentMask = false(size(grid));
                for k = 1:length(metaRow)
                    neighbors = bsxfun(@plus, neighborhoodOffsets, [metaRow(k), metaCol(k)]);
                    validNeighbors = neighbors(all(neighbors > 0 & neighbors <= gridSize, 2), :);
                    for m = 1:size(validNeighbors, 1)
                        ni = validNeighbors(m, 1);
                        nj = validNeighbors(m, 2);
                        if grid(ni, nj) == 1 && rand <= metastasisChance
                            adjacentMask(ni, nj) = true;
                        end
                    end
                end
                [adjacentRows, adjacentCols] = find(adjacentMask);
                if ~isempty(adjacentRows)
                    chosenIdx = randi(length(adjacentRows));
                    newGrid(adjacentRows(chosenIdx), adjacentCols(chosenIdx)) = 2;
                end
            end
        end

        % Update the grid
        grid = newGrid;

        % Update tumor counts
        normalTumorCount(iteration) = sum(grid(:) == 2);
        mutatedTumorCount(iteration) = sum(grid(:) == 6);
        
        % Update net survival percentage
        totalTumorCells = normalTumorCount(iteration) + mutatedTumorCount(iteration);
        netSurvival(iteration) = calculateSurvival(totalTumorCells);
        
        overviewText = sprintf(['Simulation Parameters:\n', ...
        'Break Threshold: %f\n', ...
        'Tumor Chance: %f\n', ...
        'Death Chance: %f\n', ...
        'Mutation Chance: %f\n', ...
        'Metastasis Chance: %f\n\n', ...
        'Iteration: %d\n', ...
        'Normal Tumor Cells: %d\n', ...
        'Mutated Tumor Cells: %d\n', ...
        'Net Survival: %.2f%%\n'], breakThreshold, tumorChance, deathChance, mutationChance, metastasisChance, iteration, normalTumorCount(iteration), mutatedTumorCount(iteration), netSurvival(iteration));
        set(hText, 'String', overviewText);

        if iteration == 150
            % Plot density heatmap
            figure;
            set(gcf, 'Color', 'w'); % Set background to white
            sigma = 5; % Standard deviation for the Gaussian kernel
            smoothedDensityMap = imgaussfilt(densityMap, sigma);
            
            % Normalize the density values to [0, 1]
            normalizedDensityMap = smoothedDensityMap;
            minVal = min(normalizedDensityMap(:));
            maxVal = max(normalizedDensityMap(:));
            normalizedDensityMap = (normalizedDensityMap - minVal) / (maxVal - minVal);

            % Apply custom colormap
            imagesc(normalizedDensityMap);
            % Define the custom colormap
            customCmap = [
                1 1 1;         % White for the bottom 2%
                0.22 0.21 0.87; % Blue
                0.71 0.93 0.53; % Green
                0.69 0.27 0.17  % Red
                ];
            % Number of colors
            nColors = size(customCmap, 1);
            % Interpolate colormap
            cmapInterp = interp1(linspace(0, 1, nColors), customCmap, linspace(0, 1, 256));
            % Apply the colormap
            colormap(cmapInterp);
            colorbar;
            % Set the color axis limits
            caxis([0, 1]);
            % Define custom tick marks and labels
            colorTicks = [0.02, 0.34, 0.66, 1]; % 2%, 34%, 66%, and 100%
            colorbar('Ticks', colorTicks, 'TickLabels', {'2%', '34%', '66%', '100%'});
            % Add title and axis labels
            xlabel('X', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Y', 'FontSize', 12, 'FontWeight', 'bold');
            title('Tumor Cell Density Heatmap for Day 0', 'FontSize', 14, 'FontWeight', 'bold');

            pause;
        end
        if iteration == 240
            % Plot density heatmap
            figure;
            set(gcf, 'Color', 'w'); % Set background to white
            sigma = 5; % Standard deviation for the Gaussian kernel
            smoothedDensityMap = imgaussfilt(densityMap, sigma);
            
            % Normalize the density values to [0, 1]
            normalizedDensityMap = smoothedDensityMap;
            minVal = min(normalizedDensityMap(:));
            maxVal = max(normalizedDensityMap(:));
            normalizedDensityMap = (normalizedDensityMap - minVal) / (maxVal - minVal);

            % Apply custom colormap
            imagesc(normalizedDensityMap);
            % Define the custom colormap
            customCmap = [
                1 1 1;         % White for the bottom 2%
                0.22 0.21 0.87; % Blue
                0.71 0.93 0.53; % Green
                0.69 0.27 0.17  % Red
                ];
            % Number of colors
            nColors = size(customCmap, 1);
            % Interpolate colormap
            cmapInterp = interp1(linspace(0, 1, nColors), customCmap, linspace(0, 1, 256));
            % Apply the colormap
            colormap(cmapInterp);
            colorbar;
            % Set the color axis limits
            caxis([0, 1]);
            % Define custom tick marks and labels
            colorTicks = [0.02, 0.34, 0.66, 1]; % 2%, 34%, 66%, and 100%
            colorbar('Ticks', colorTicks, 'TickLabels', {'2%', '34%', '66%', '100%'});
            % Add title and axis labels
            xlabel('X', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Y', 'FontSize', 12, 'FontWeight', 'bold');
            title('Tumor Cell Density Heatmap for Day 3', 'FontSize', 14, 'FontWeight', 'bold');

            pause;
        end
    end

    % Save tumor counts to CSV
    tumorData = table((1:numIterations)', normalTumorCount', mutatedTumorCount', 'VariableNames', {'Iteration', 'NormalTumorCount', 'MutatedTumorCount'});
    writetable(tumorData, 'tumor_counts.csv');

    % Plot tumor growth over time
    figure;
    set(gcf, 'Color', 'w'); % Set background to white
    hold on; % Hold on to add multiple lines to the same plot
    plot(1:numIterations, normalTumorCount, 'b-', 'LineWidth', 2);
    plot(1:numIterations, mutatedTumorCount, 'r-', 'LineWidth', 2);
    hold off; % Release the hold
    legend('Normal Tumor Cells', 'Mutated Tumor Cells');
    xlabel('Iteration', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Number of Cells', 'FontSize', 12, 'FontWeight', 'bold');
    title('Tumor Growth Over Time', 'FontSize', 14, 'FontWeight', 'bold');

    % Plot net survival over time with error range
    figure;
    set(gcf, 'Color', 'w'); % Set background to white
    % Calculate the error range (example: 10% of the tumor cell count)
    errorRange = 0.001 * (normalTumorCount + mutatedTumorCount);
    % Plot the survival curve
    plot(1:numIterations, netSurvival, 'r-', 'LineWidth', 2);
    hold on;
    % Plot the error range
    fill([1:numIterations, fliplr(1:numIterations)], ...
        [netSurvival - errorRange, fliplr(netSurvival + errorRange)], ...
        'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    % Set y-axis limits
    ylim([0, 100]);
    xlabel('Iteration', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Net Survival Percentage (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Net Survival Percentage Over Time', 'FontSize', 14, 'FontWeight', 'bold');
    hold off;
        
    % Plot final density heatmap
    figure;
    set(gcf, 'Color', 'w'); % Set background to white
    sigma = 5; % Standard deviation for the Gaussian kernel
    smoothedDensityMap = imgaussfilt(densityMap, sigma);
        
    % Normalize the density values to [0, 1]
    normalizedDensityMap = smoothedDensityMap;
    minVal = min(normalizedDensityMap(:));
    maxVal = max(normalizedDensityMap(:));
    normalizedDensityMap = (normalizedDensityMap - minVal) / (maxVal - minVal);

    % Apply custom colormap
    imagesc(normalizedDensityMap);
    % Define the custom colormap
    customCmap = [
        1 1 1;         % White for the bottom 2%
        0.22 0.21 0.87; % Blue
        0.71 0.93 0.53; % Green
        0.69 0.27 0.17  % Red
        ];
    % Number of colors
    nColors = size(customCmap, 1);
    % Interpolate colormap
    cmapInterp = interp1(linspace(0, 1, nColors), customCmap, linspace(0, 1, 256));
    % Apply the colormap
    colormap(cmapInterp);
    colorbar;
    % Set the color axis limits
    caxis([0, 1]);
    % Define custom tick marks and labels
    colorTicks = [0.02, 0.34, 0.66, 1]; % 2%, 34%, 66%, and 100%
    colorbar('Ticks', colorTicks, 'TickLabels', {'2%', '34%', '66%', '100%'});
    % Add title and axis labels
    xlabel('X', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize', 12, 'FontWeight', 'bold');
    title('Tumor Cell Density Heatmap for Day 5', 'FontSize', 14, 'FontWeight', 'bold');
end

function netSurvival = calculateSurvival(totalTumorCells)
    % Parameters for the survival model
    initialSurvivalRate = 100; % 100% initial survival rate
    maxTumorCells = 20000; % Maximum number of tumor cells where survival rate drops to 0%

    % Calculate net survival rate
    netSurvival = initialSurvivalRate * (1 - totalTumorCells / maxTumorCells);

    % Ensure survival does not go below 0%
    netSurvival(netSurvival < 0) = 0;
end
