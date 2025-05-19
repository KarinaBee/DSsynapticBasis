load pTable.mat

% Iterate over all recordings
for i = 1:height(pTable)
    
    % Extract peakCurrCorr and VhCorr values for this recording
    I = pTable.peakCurr{i};  % Current values (in pA)
    V = [-80; -70;-60; -50; -40];
    % VhCorr{i};       
    % I = pTable.peakCurrCorr{i};  % Current values (in pA)
    % V = pTable.VhCorr{i};        % Voltage values (in mV)

    % Ensure the data isn't empty or NaN-filled before proceeding
    if isempty(I) || isempty(V) || all(isnan(I)) || all(isnan(V))
        warning('Skipping row %d due to empty or NaN values.', i);
        continue;
    end

    % Ensure both I and V are column vectors
    I = I(:);  % Convert to column vector
    V = V(:);  % Convert to column vector

    % Use up to the first 6 data points, or as many as are available
    numPoints = min(5, length(I));  % Determine how many points are available
    I = I(1:numPoints);  % Take up to the first 6 current values
    V = V(1:numPoints);  % Take up to the first 6 voltage values

    % Perform linear regression using fitlm
    mdl = fitlm(V, I);  
    R_squared = mdl.Rsquared.Ordinary;  % R-squared value
    y_est = predict(mdl, V);  % Predicted fitted values

    % Extract slope (conductance) and intercept from the model
    slope = mdl.Coefficients.Estimate(2);  % Slope (Conductance, nS)
    intercept = mdl.Coefficients.Estimate(1);  % Intercept

    % Save the regression results into pTable
    pTable.P1(i) = slope;
    pTable.P2(i) = intercept;
    pTable.RSquared(i) = R_squared;
    pTable.gGaba(i) = slope;  % Store conductance in nS

    % Plot data points and the fitted line
    figure;
    plot(V, I, 'o', 'Color', 'b', 'DisplayName', 'Data Points');
    hold on;
    plot(V, y_est, '--', 'Color', 'r', 'LineWidth', 2);

    % Add labels and title
    xlabel('Voltage (mV)');
    ylabel('Current (pA)');
    % ylim([-200 200]);

    % Get experiment date and acquisition number for title
    expDate = pTable.expDate(i);
    acqNum = pTable.aquisitionNumber(i);
    sacSide = pTable.SACside(i);
    geno = pTable.genotype(i);
    gfp = pTable.GFPLabel(i);
    title([expDate geno gfp acqNum sacSide]);
    legend(['R^2 = ', num2str(R_squared, '%.4f'), ...
            ', gGABA = ', num2str(slope, '%.3f'), ' nS'], 'Location', 'best');
end

% Save the updated pTable to a file
% cd 'C:\Users\Karina Bistrong\Dropbox\Mac\Desktop\DATA\Ephys\1_compiledResults\Pairs'
% save('updated_pTable.mat', 'pTable'); 