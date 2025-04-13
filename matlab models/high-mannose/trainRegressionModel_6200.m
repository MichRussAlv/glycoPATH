% Define the base name and suffix pattern
varBase = 'filtered_HM_genes_transposed_test_';
test = 'H7_N2_F0_S0';

% Construct the variable name
varName = [varBase, test];

% Access the variable from the workspace
if evalin('base', ['exist(''', varName, ''', ''var'')'])
    % Retrieve the variable value from the base workspace
    data = evalin('base', varName);
else
    error('Variable %s does not exist in the workspace.', varName);
end

% Make predictions using the trained model,update glycan model
yfit = trainedModel7200.predictFcn(data);

% Define the new column name for predictions
predColumnName = 'Predictions';

% Add predictions as a new column in the table
data.(predColumnName) = yfit;

% Define the filename for the CSV output
csvFilename = [varBase, test, '_predictions.csv'];

% Export the updated table to a CSV file
writetable(data, csvFilename);
