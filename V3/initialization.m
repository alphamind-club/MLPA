close all

macroParam = [0.001, 0.05, 0.95, 0.01, 0.000005];

% Step 1: Read the first column of the CSV file for gene pathway names
filename = 'pathways.csv'; % Adjust the path as necessary
opts = detectImportOptions(filename, 'NumHeaderLines', 0, 'VariableNamingRule', 'preserve'); % Adjust based on your CSV file's header
opts.SelectedVariableNames = [1]; % Assuming the gene pathways are in the first column
genePathwaysTable = readtable(filename, opts);

% Convert table to cell array of strings for gene names
genePathwayNames = table2cell(genePathwaysTable);

% Initialize an empty cell array to hold interaction lists
interactionLists = cell(size(genePathwayNames, 1), 1);
for i = 1:size(genePathwayNames, 1)
    interactionLists{i} = genePathwayNames(i, :); % Assign the row as an interaction list
end
interactionLists = generateInteractionListsFromCSV();

% Step 2 and 3: Create a resizable GUI dialog for each pathway and collect inputs
microParam = createResizableGUIForGenePathways(genePathwayNames);
disp('Micro Parameters:');
disp(microParam);

for i = 1:length(microParam)
    if microParam(i) == 1  % Check if the gene pathway is expressed
        for j = 1:length(macroParam)  % Iterate through each macro parameter
            interaction = interactionLists{i}{j};  % Get the interaction for this gene pathway and macro parameter
            switch interaction
                case '+'  % Activation interaction
                    macroParam(j) = macroParam(j) + 10;
                case '-'  % Inhibition interaction
                    macroParam(j) = macroParam(j) - 10;
                case '?'  % Unknown interaction, consider custom logic or leave as is
                    % Implement any custom logic for '?' if needed
                case 'NA'  % No action needed
                    % Do nothing for 'NA'
            end
        end
    end
end
disp('Macro Parameters:');
disp(macroParam);

camodel(macroParam(1), macroParam(2), macroParam(3), macroParam(4), macroParam(5));

function interactionLists = generateInteractionListsFromCSV()
    % Step 1: Read CSV File Excluding the First Column
    filename = 'pathways.csv'; % Specify your CSV file name
    opts = detectImportOptions(filename, 'VariableNamingRule', 'preserve'); % Automatically detect import options
    opts.SelectedVariableNames = opts.VariableNames(2:end); % Exclude the first column
    interactionData = readtable(filename, opts); % Read the CSV file excluding the first column

    % Convert table to cell array for easier processing
    interactionArray = table2cell(interactionData);

    % Step 2: Construct Interaction Lists
    % Initialize an empty cell array to hold interaction lists
    interactionLists = cell(size(interactionArray, 1), 1);

    % Iterate through each row of interactionArray to create interaction lists
    for i = 1:size(interactionArray, 1)
        interactionLists{i} = interactionArray(i, :); % Assign the row as an interaction list
    end
end

function microParam = createResizableGUIForGenePathways(genePathwayNames)
    microParam = []; % Initialize to empty; will be updated by the submit button callback
    figHeight = max(550, 50 + 25 * length(genePathwayNames)); 
    fig = figure('Name', 'Gene Pathway Expression', ...
                 'Position', [100, 100, 400, figHeight], ...
                 'Resize', 'on', ...
                 'CloseRequestFcn', @closeFig); 

    numPathways = length(genePathwayNames);
    startY = figHeight - 40; 
    checkboxes = gobjects(numPathways, 1); 
    for i = 1:numPathways
        checkboxes(i) = uicontrol('Style', 'checkbox', ...
                                  'Parent', fig, ...
                                  'String', genePathwayNames{i}, ...
                                  'Value', 0, ...
                                  'Position', [20, startY - (i-1) * 25, 360, 20], ...
                                  'HandleVisibility', 'off');
    end

    submitButton = uicontrol('Style', 'pushbutton', ...
                             'String', 'Submit', ...
                             'Position', [(400-200)/2, 15, 200, 30], ...
                             'Callback', @submitButtonCallback); 

    uiwait(fig); % Wait for GUI interaction

    function submitButtonCallback(~, ~)
        microParamTemp = zeros(1, numPathways); 
        for j = 1:numPathways
            microParamTemp(j) = get(checkboxes(j), 'Value'); 
        end
        microParam = microParamTemp; % Update the outer variable
        uiresume(fig); % Allow the script to continue
        close(fig); 
    end

    function closeFig(~, ~)
        uiresume(fig); 
        delete(fig);
    end
end
