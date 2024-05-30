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
    minTumorCellCountForMetastasis = 300;
    treatment_time = [100, 200, 300];

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
    startCol = randi([1, gridSize - vesselWidth + 1]);
    vesselMask = false(gridSize, gridSize);
    vesselMask(startRow:(startRow + vesselHeight - 1), startCol:(startCol + vesselWidth - 1)) = true;
    grid(vesselMask) = 5; % Assign a new value (e.g., 5) for blood vessel

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
    growthRateMap(polarizedMap > 0.6) = 0.02;      % Values above 0.6 set to 1
    growthRateMap(polarizedMap > 0.3 & polarizedMap <= 0.6) = 0.02; % Values between 0.3 and 0.6 set to 0.5
    growthRateMap(polarizedMap <= 0.3) = 1;     % Values below or equal to 0.3 set to 0

    % Display the growth rate map
    figure('Name', 'Polarized Map', 'NumberTitle', 'off');
    imagesc(polarizedMap);
    colorbar;
    title('Polarized Map');
    xlabel('X');
    ylabel('Y');

    figure('Name', 'Growth bdRate Map', 'NumberTitle', 'off');
    imagesc(growthRateMap);
    colorbar;
    title('Growth Rate Map');
    xlabel('X');
    ylabel('Y');

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
    numIterations = 500;
    numberOfBreaks = 0;
    normalTumorCount = zeros(1, numIterations);
    mutatedTumorCount = zeros(1, numIterations);
    maxDistance = abs(center(2) - startCol);
    metastasis = 0;

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
        title(['Iteration: ' num2str(iteration)]);
        drawnow;

        % Update grid
        newGrid = grid;

        % Tumor growth
        growMask = false(size(grid));
        for i = 1:gridSize
            for j = 1:gridSize
                if grid(i, j) == 2 || grid(i, j) == 6
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
                
                %create a break
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

                %select the closest tumor cell to grow to
                tumorDistance = 1e+30;
                if grid(i,j) == 4 && tumorCount == 0
                    for a = 1:gridSize
                        for b =1:gridSize
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

                %grow towards the closest tumor cell
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

                    %make sure it doesnt grow outside the box
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

                %stop if touch tumor
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
        overviewText = sprintf(['Simulation Parameters:\n', ...
        'Break Threshold: %f\n', ...
        'Tumor Chance: %f\n', ...
        'Death Chance: %f\n', ...
        'Mutation Chance: %f\n', ...
        'Metastasis Chance: %f\n\n', ...
        'Iteration: %d\n', ...
        'Normal Tumor Cells: %d\n', ...
        'Mutated Tumor Cells: %d\n'], breakThreshold, tumorChance, deathChance, mutationChance, metastasisChance, iteration, normalTumorCount(iteration), mutatedTumorCount(iteration));
    set(hText, 'String', overviewText);
    end

    % Plot tumor growth over time
    figure;
    plot(1:numIterations, normalTumorCount, 'b-', 1:numIterations, mutatedTumorCount, 'r-');
    legend('Normal Tumor Cells', 'Mutated Tumor Cells');
    xlabel('Iteration');
    ylabel('Number of Cells');
    title('Tumor Growth Over Time');
end
