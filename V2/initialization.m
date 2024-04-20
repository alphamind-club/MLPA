% Define Micro Parameters, I(i)
microParams = [10, 5, 15, 20]; % Example values for I(1), I(2), I(3), and I(4)

% Define Macro Parameters, A(j)
breakThreshold = 0.001; % Example value for breakThreshold
tumorChance = 0.02; % Example value for tumorChance
deathChance = 0.95; % Example value for deathChance
mutationChance = 0.01; % Example value for mutationChance
metastasisChance = 0.0005; % Example value for metastasisChance

macroParams = [breakThreshold, tumorChance, deathChance, mutationChance, metastasisChance];

% Define interactions from micro to macro parameters
interactions = {
    {'NA', 'NA', '+', 'NA'};  % For Macro 1: Only Micro 3 activates
    {'-', 'NA', 'NA', '+'};   % For Macro 2: Micro 1 inhibits, Micro 4 activates
    {'NA', '?', 'NA', 'NA'};   % For Macro 3: Unknown association with Micro 2
    {};
    {}
};

% Now, let's run the function and print the output
outputs = MicroMacroInteraction(microParams, macroParams, interactions);

% Display the results
disp('Output Parameters:');
disp(outputs);

camodel(macroParams(1), macroParams(2), macroParams(3), macroParams(4), macroParams(5))

function outputs = MicroMacroInteraction(microParams, macroParams, interactions)
    % microParams: Array of micro level parameters I
    % macroParams: Array of macro level parameters A
    % interactions: Cell array of interaction types between micro and macro levels

    % Initialize output parameters array
    outputs = zeros(1, length(macroParams));

    % Loop through each macro parameter
    for j = 1:length(macroParams)
        % Reset influence for each macro parameter
        influence = 0;

        for i = 1:length(interactions{j})
            % Get the interaction type
            interactionType = interactions{j}{i};
 
            % Calculate influence based on interaction type
            switch interactionType
                case '+'
                    influence = influence + microParams(i); % Activation
                case '-'
                    influence = influence - microParams(i); % Inhibition
                case '?'
                    % Unknown association, could be a custom implementation
                    % For simplicity, we'll leave it as is or you can define a rule
                case 'NA'
                    % No association, do nothing
                otherwise
                    % Handle unexpected interaction type
                   fprintf('Unknown interaction type: %s\n', interactionType);
            end
        end

        % Assign the influenced value to the output parameter
        if (influence ~= 0)
            outputs(j) = macroParams(j) * influence;
        else
            outputs(j) = macroParams(j);
        end
        disp(outputs(j))
    end
end
