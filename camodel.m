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

    % Visualization setup
    figure;
    colormap([1 1 1; 0 0 1; 1 0.75 0.7; 0.5 0 0; 1 0.75 0.7; 0.5 0.5 1]); % Color map

    % Simulation parameters
    numIterations = 500;
    numberOfBreaks = 0;
    normalTumorCount = zeros(1, numIterations);
    mutatedTumorCount = zeros(1, numIterations);
    maxDistance = abs(center(2) - startCol);
    metastasis = 0;

    % Run simulation
    for iteration = 1:numIterations
        % Visualization
        image(grid);
        title(['Iteration: ' num2str(iteration)]);
        drawnow;

        %update grid
        newGrid = grid;
        for i = 1:gridSize
            for j = 1:gridSize
                neighborhood = grid(max(1, i-neighborhoodSize):min(gridSize, i+neighborhoodSize), max(1, j-neighborhoodSize):min(gridSize, j+neighborhoodSize));
                emptySpace = sum(neighborhood(:) == 1);
                tumorCount = sum(neighborhood(:) == 2) + sum(neighborhood(:) == 6);
                vesselCount = sum(neighborhood(:) == 3);

                growhood = grid(max(1, i-growhoodSize):min(gridSize, i+growhoodSize), max(1, j-growhoodSize):min(gridSize, j+growhoodSize));
                canGrow = sum(growhood(:) == 3);

                %mutate
                if grid(i, j) == 2 && rand <= mutationChance
                    newGrid(i, j) = 6;
                end

                %grow tumor
                if grid(i, j) == 1 && canGrow > 0 && tumorCount > 0 && rand <= tumorChance
                    newGrid(i, j) = 2;
                elseif grid(i, j) == 1 && tumorCount > 0 && rand <= tumorChance / 5
                    newGrid(i, j) = 2;
                end

                if grid(i, j) == 1 && canGrow > 0 && tumorCount > 0 && rand <= (tumorChance / 1.5)
                    newGrid(i, j) = 2;
                elseif grid(i, j) == 1 && tumorCount > 0 && rand <= ((tumorChance / 5) / 1.5)
                    newGrid(i, j) = 2;
                end

                %drug administration
                if ismember(iteration, treatment_time)
                    if grid(i, j) == 2 && rand <= deathChance
                        newGrid(i, j) = 1;
                    elseif grid(i, j) == 6 && rand <= (deathChance / 2)
                        newGrid(i, j) = 1;
                    end
                end

                %create a break
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

                %stop if touches tumor
                if grid(i,j) == 4 && tumorCount ~= 0
                    newGrid(i, j) = 3;
                    metastasis = 1;
                end
            end

            %metasis
            metastasisPoint = [randi([startRow, vesselHeight]), 0];
            if rand >= 0.5
                metastasisPoint(2) = startCol - 1;
            else
                metastasisPoint(2) = startCol + vesselWidth;
            end
            if metastasis == 1 && rand <= metastasisChance && sum(grid(:) == 2) + sum(grid(:) == 6) >= minTumorCellCountForMetastasis && grid(metastasisPoint(1), metastasisPoint(2)) == 1
                newGrid(metastasisPoint(1), metastasisPoint(2)) = 2;
            end

            normalTumorCount(iteration) = sum(grid(:) == 2);
            mutatedTumorCount(iteration) = sum(grid(:) == 6);
        end
        grid = newGrid;
    end

    % Plot the tumor cell counts over time
    figure;
    plot(1:numIterations, normalTumorCount, 'b-', 'LineWidth', 2, 'DisplayName', 'Normal Tumor Cells');
    hold on;
    plot(1:numIterations, mutatedTumorCount, 'r-', 'LineWidth', 2, 'DisplayName', 'Mutated Tumor Cells');
    xlabel('Iteration');
    ylabel('Cell Count');
    title('Normal vs Mutated Tumor Cell Count Over Time');
    legend('show');
end