% this code takes pTable, with the updated values from running IVfit4, and
% extracts each conductance value 
% Then it
% extracts the same list of the SAC position, and separates each conductance
% with the corresponding SAC side - in order to plot all the conductances
% for both null and pref side recordings 

% Load the table with updated conductance values
load pTableMaster.mat 

% RIGHT NOW ITS RAW
% pTable = pTable(strcmp(pTable.genotype,"WT") & strcmp(pTable.GFPLabel,"DRD4") ,:);
% pTable = pTable(strcmp(pTable.genotype,"B2") & strcmp(pTable.GFPLabel,"DRD4") ,:);
pTable = pTable(strcmp(pTable.genotype,"WT") & strcmp(pTable.GFPLabel,"HB9") ,:);
% pTable = pTable(strcmp(pTable.genotype,"B2") & strcmp(pTable.GFPLabel,"HB9") ,:);

% Extract relevant data: Conductance (gGABA), SAC side, and distance
pData = pTable.gGaba;      % Conductance values
pSide = pTable.SACside;    % SAC side labels ('pref' or 'null')
pDist = pTable.distance;   % Distance values
pOverlap = pTable.overlap; % Overlap values

% Ensure data are column vectors
pData = pData(:);
pSide = pSide(:);
pDist = pDist(:);
pOverlap = pOverlap(:);

% Separate gGABA and distance values by SAC side (Pref/Null)
pref_gGaba = pData(strcmp(pSide, 'pref'));
pref_Dist = pDist(strcmp(pSide, 'pref'));
pref_Dist = str2double(pref_Dist);

null_gGaba = pData(strcmp(pSide, 'null'));
null_Dist = pDist(strcmp(pSide, 'null'));
null_Dist = str2double(null_Dist);

% Plot the conductance values for both 'pref' and 'null' sides
figure('Name', 'GABA Conductance');
tiledlayout(1, 2);  % Create two side-by-side plots
% Plot the 'pref' side conductance values
nexttile;
plotSpread(pref_gGaba, 'showMM', 5);  % Spread plot with mean/median markers
xlabel('Pref');
ylabel('Conductance (nS)');
ylim([0 10]);
title('Pref Side GABA Conductance');

% Plot the 'null' side conductance values
nexttile;
plotSpread(null_gGaba, 'showMM', 5);  % Spread plot with mean/median markers
xlabel('Null');
ylabel('Conductance (nS)');
ylim([0 10]);
title('Null Side GABA Conductance');

% --- Linear Regression on Distance ---
% Perform linear regression for pref and null sides for Distance
mdl_pref_dist = fitlm(pref_Dist, pref_gGaba);
mdl_null_dist = fitlm(null_Dist, null_gGaba);

% --- Linear Regression on Overlap ---
% Filter valid rows where overlap values exist (not NaN)
validIdx = ~isnan(pOverlap);
valid_gGaba = pData(validIdx);
valid_Side = pSide(validIdx);
valid_Overlap = pOverlap(validIdx);


% Separate by SAC side (Pref/Null)
pref_Overlap = valid_Overlap(strcmp(valid_Side, 'pref'));
pref_gGaba_overlap = valid_gGaba(strcmp(valid_Side, 'pref'));

null_Overlap = valid_Overlap(strcmp(valid_Side, 'null'));
null_gGaba_overlap = valid_gGaba(strcmp(valid_Side, 'null'));


%% Plot 2: gGABA as a Function of Distance

figure('Name', 'gGABA vs Distance');

% Plot for Pref side
plot(pref_Dist, pref_gGaba, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Pref');
hold on;
% Plot the fitted regression line for Pref side
plot(pref_Dist, mdl_pref_dist.Fitted, 'b-', 'LineWidth', 2, 'DisplayName', 'Pref Fit');

% Add R^2 value for Pref side
text_x_pref = min(pref_Dist) + 0.1 * range(pref_Dist); % Position on x-axis
text_y_pref = max(pref_gGaba) - 0.1 * range(pref_gGaba); % Position on y-axis
text(text_x_pref, text_y_pref, sprintf('R^2 = %.2f', mdl_pref_dist.Rsquared.Ordinary), 'Color', 'b', 'FontSize', 10);

% Plot for Null side
plot(null_Dist, null_gGaba, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Null');
% Plot the fitted regression line for Null side
plot(null_Dist, mdl_null_dist.Fitted, 'r-', 'LineWidth', 2, 'DisplayName', 'Null Fit');

% Add R^2 value for Null side
text_x_null = min(null_Dist) + 0.1 * range(null_Dist); % Position on x-axis
text_y_null = max(null_gGaba) - 0.1 * range(null_gGaba); % Position on y-axis
text(text_x_null, text_y_null, sprintf('R^2 = %.2f', mdl_null_dist.Rsquared.Ordinary), 'Color', 'r', 'FontSize', 10);


% Add labels, legend, and title
xlabel('Distance (\mum)');
ylabel('Conductance (nS)');
title('gGABA as a Function of Distance');
legend('Location', 'best');
xlim([10 130]);
ylim([0 10]);
grid off;
hold off;

% Plot conductance vs overlap
figure('Name', 'gGABA vs Overlap');

% Plot for Pref side
plot(pref_Overlap, pref_gGaba_overlap, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Pref');
hold on;
% Perform linear regression on Pref side for overlap
mdl_pref_overlap = fitlm(pref_Overlap, pref_gGaba_overlap);
% Plot the fitted regression line for Pref side
plot(pref_Overlap, mdl_pref_overlap.Fitted, 'b-', 'LineWidth', 2, 'DisplayName', 'Pref Fit');

% Add R^2 value for Pref side
text_x_pref_overlap = min(pref_Overlap) + 0.1 * range(pref_Overlap); % Position on x-axis
text_y_pref_overlap = max(pref_gGaba_overlap) - 0.1 * range(pref_gGaba_overlap); % Position on y-axis
text(text_x_pref_overlap, text_y_pref_overlap, sprintf('R^2 = %.2f', mdl_pref_overlap.Rsquared.Ordinary), 'Color', 'b', 'FontSize', 10);

% Plot for Null side
plot(null_Overlap, null_gGaba_overlap, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Null');
% Perform linear regression on Null side for overlap
mdl_null_overlap = fitlm(null_Overlap, null_gGaba_overlap);
% Plot the fitted regression line for Null side


% Plot the fitted regression line for Null side (spanning min to max overlap)
x_null_overlap_fit = linspace(min(null_Overlap), max(null_Overlap), 100);  % Generate x values
y_null_overlap_fit = predict(mdl_null_overlap, x_null_overlap_fit');  % Compute corresponding y values
plot(x_null_overlap_fit, y_null_overlap_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Null Fit');
% plot(null_Overlap, mdl_null_overlap.Fitted, 'r-', 'LineWidth', 2, 'DisplayName', 'Null Fit');

% Add R^2 value for Null side
text_x_null_overlap = min(null_Overlap) + 0.1 * range(null_Overlap); % Position on x-axis
text_y_null_overlap = max(null_gGaba_overlap) - 0.1 * range(null_gGaba_overlap); % Position on y-axis
text(text_x_null_overlap, text_y_null_overlap, sprintf('R^2 = %.2f', mdl_null_overlap.Rsquared.Ordinary), 'Color', 'r', 'FontSize', 10);

% Add labels, legend, and title
xlabel('Overlap');
ylabel('Conductance (nS)');
title('gGABA vs Overlap (Pref and Null)');
legend('Location', 'best');
xlim([1000 17000]);
ylim([0 5]);
grid off;
hold off;
