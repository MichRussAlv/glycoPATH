% Assuming trainedModel is your trained model and validationData is your validation dataset
% Make predictions on validation data
yPred = trainedModel.predictFcn(validationData);

% True labels for validation data
trueLabels = validationData.TrueLabels;  % Adjust if needed

% Calculate RMSE
rmse = sqrt(mean((yPred - trueLabels).^2));

% Create a table for results
resultsTable = table(yPred, trueLabels, 'VariableNames', {'Predictions', 'TrueLabels'});

% Add RMSE to the table
resultsTable.RMSE = repmat(rmse, height(resultsTable), 1);  % Add RMSE as a column

% Define the filename for the CSV output
csvFilename = 'validation_results_with_rmse.csv';

% Export the table to a CSV file
writetable(resultsTable, csvFilename);
