%% Build the neuron table for PAIRS

%This code just loads the data into pTable WITH PLOTS. See commented sections below 


% This code needs an excel spreadsheet with metadata in it.
% Right now, it will read 'PairsExperimentSheet.xlsx'
% Every row in that excel spreadsheet is a single recording
% The buildTable code will run all the recording
%
% The end result of this code is 'pTable.mat'
% This is the table where each row is a recording
% Other codes will load this table and run analysis on it
% no functions for this code

% in order to get the conductance value you must run IVfit4 --> SAVE -->
% run plotgGABA --> SAVE

%% Load the table
data_guide_name = '';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'ExperimentDate','macPath','path','FileName','Sex','Age', 'Genotype', 'GFPLabel',...
    'Eye', 'Quadrant','SACside','StimSide','StimDelay','CellNumb','Distance'},'string');

expTable = readtable(data_guide_name, opts');

[num_files, dummy] = size(expTable); %determines number of rows (iterations for loop will go through)

%% Intialize the neuron table (called pTable)

% Initialize the neuron table (called pTable)
totalRecs = height(expTable);

pTable = table('Size', [totalRecs 18], 'VariableTypes', ...
    {'double','double','string','string','string','string','string','string','string','string', ...
    'string','string','double','double', 'cell', 'cell', 'cell', 'cell'});

pTable.Properties.VariableNames = {'expDate', 'aquisitionNumber', 'sex', 'age', 'genotype', ...
    'GFPLabel', 'eye', 'quad', 'SACside', 'stimSide', 'cellNumb', 'distance', 'overlap','holdVol', ...
    'VhCorr', 'peakCurr', 'peakCurrCorr', 'C2corrected'};

pTable.C1norm = NaN(totalRecs, 20000); % Optional: Keep this as-is for fixed size

neuronCounter = 1; % Initialize the neuron counter

for i = 1:num_files
   % Loop through all recordings
    % Change directory to the folder where data is stored for this recording
    % cd(expTable.path(i));

    cd(expTable.macPath(i));
    % Load the ABF file data
    [d, si] = abfload(expTable.FileName(i));
    C1_data = squeeze(d(:, 1, :)); % Presynaptic data
    C2_data = squeeze(d(:, 2, :)); % Postsynaptic data

    % Retrieve metadata for this recording
    holdVol = expTable.holdVol(i);
    if holdVol == 9
        Vhold = [-80; -70; -60; -50; -40; -30; -20; -10; 0];
    elseif holdVol == 6
        Vhold = [-80; -70; -60; -50; -40; -30];
    elseif holdVol == 5
        Vhold = [-80; -70; -60; -50; -40];
    end

    % Initialize arrays for storing results for this row
    peakCurrList = NaN(1, holdVol);
    peakCurrCorrList = NaN(1, holdVol);
    VhCorrList = NaN(1, holdVol);
    C2correctedList = cell(1, holdVol); % Store corrected C2 data for each hold voltage
    RsList = NaN(1,holdVol);
    RinList = NaN(1, holdVol);
    totalRList = NaN(1,holdVol); 
    C1norm_all = []; % This will store all C1norm traces for averaging

    showPlots = true; % Set to false to suppress plots
    % showPlots = false; % Set to false to suppress plots

   if showPlots    % Create a new figure for plotting the current traces
        figure;
        hold on;
        legendEntries = cell(1, holdVol); % Preallocate legend entries
   end 
    % Set the appropriate baselines and time windows based on StimDelay
    if expTable.StimDelay(i) == "150"
        DSbaseline = 1800:1850;
        sacBaseLine = 500:1500; 
        timeOfPeak = 1800:2000;
        stimWindow = 1500:3800;
    elseif expTable.StimDelay(i) == "200"
        DSbaseline = 2300:2350;
        sacBaseLine = 1000:2000; 
        timeOfPeak = 2400:2800;
        stimWindow = 2200:4200;
    elseif expTable.StimDelay(i) == "500"
        DSbaseline = 2800:2850;
        sacBaseLine = 1000:2000; 
        timeOfPeak = 3200:4000;
        stimWindow = 3000:5000;
    end

        % Loop over each holding voltage
    for j = 1:min(5, holdVol) % Only iterate through the first 6 hold voltages
        % Normalize baseline currents for C2 and C1
        C2 = C2_data(:, j);  
        C2norm = C2 - mean(C2(DSbaseline)); % Subtract baseline
        C2flip = C2norm'; % Transpose for plotting (optional)

        C1 = C1_data(:, j);  
        C1norm = C1 - mean(C1(sacBaseLine)); % Subtract SAC baseline
        C1norm = C1norm'; 

         % Store C1norm data for averaging
        C1norm_all(j, :) = C1norm(stimWindow); % Store only the stimWindow part


        % Calculate series resistance (Rs) and input resistance (Rin)
        ssCurr = abs(mean(C2norm(14400:14800))) - abs(mean(C2norm(13500:14000)));
        ssCurr = ssCurr * 10^-12; % pA to A
        dV = 5 / 1000; % mV to V
        Rin = abs(dV / ssCurr); 

        baseline_window = 13500:14000;
        peak_window = 14000:14500;

        baselineCurr = mean(C2norm(baseline_window));
        peakCurr = min(C2norm(peak_window)); % Use max if expecting outward currents
        dCurr = abs(peakCurr - baselineCurr);

        dCurr = dCurr * 10^-12; % pA to A
        Rs = abs(dV / dCurr); 

        % Correct C2 using total resistance
        totalR = (Rin + Rs) / Rin;
        C2corrected = C2norm * totalR;
        C2correctedList{j} = C2corrected(stimWindow); % Save corrected trace

        RsList(j) = Rs; 
        RinList(j) = Rin;
        totalRList(j) = totalR; 

        if showPlots
          % Plot the corrected C2 data
        plot(C2corrected(stimWindow));
        legendEntries{j} = ['Vhold = ' num2str(Vhold(j)) ' mV']; % Add legend entry
        end 

        % Find the peak current and corrected peak current
        if mean(C2corrected(timeOfPeak)) < 0
    % Inward current
    peakCurrCorr = min(C2corrected(timeOfPeak));
    peakCurr = min(C2norm(timeOfPeak));
        else
    % Outward current
    peakCurrCorr = max(C2corrected(timeOfPeak));
    peakCurr = max(C2norm(timeOfPeak));
        end

        % Calculate corrected holding voltage (VhCorr)
        mVh = Vhold(j);
        Vh = mVh / 1000; % mV to V
        peakCurrA = peakCurr * 10^-12; % pA to A
        VhCorr = Vh - (peakCurrA * Rs);
        VhCorrList(j) = VhCorr * 1000; % Convert back to mV

        % Store peak currents and corrected voltage for this voltage step
        peakCurrList(j) = peakCurr;
        peakCurrCorrList(j) = peakCurrCorr;
    end

    % Calculate the average C1norm across all holding voltages
    C1normAvg = mean(C1norm_all, 1); % Average across rows (holding voltages)

       % Get the expDate for this set
         ExpDate = expTable.ExperimentDate(i);
        geno = expTable.Genotype(i);
         gfp = expTable.GFPLabel(i);
        AcqNum = expTable.AquisitionNumber(i);
         SACside = expTable.SACside(i);
    % Update plot to display only the first 6 data points
   if showPlots
    hold on;
    xlabel('Time (ms)');
    ylabel('Current (pA)');
    legend(legendEntries(1:min(5, holdVol))); % Use only legends for the first 6 voltages
    title([ExpDate geno gfp AcqNum SACside]);
    ylim([-150 150]);
 hold off;
 %    figure;
 %    plot(Vhold(1:j), peakCurrList(1:j), '-o', 'DisplayName', 'Raw Currents');
 %    hold on;
 %    plot(VhCorrList(1:j), peakCurrCorrList(1:j), '-x', 'DisplayName', 'Corrected Currents');
 %    xlabel('Holding Voltage (mV)');
 %    ylabel('Peak Current (pA)');
 %    ylim([-200 200]);
   end
%    if showPlots
%        % Plot the averaged C1norm
% figure;
% plot(C1normAvg);
% ylim([-1000 1000]);
% xlabel('Time (ms)');
% ylabel('Averaged Normalized Current (pA)');
% title('Average C1norm Across Holding Voltages');
%    end 

    % Store all data in pTable for this recording
    pTable.expDate(neuronCounter) = expTable.ExperimentDate(i);
    pTable.aquisitionNumber(neuronCounter) = expTable.AquisitionNumber(i);
    pTable.sex(neuronCounter) = expTable.Sex(i);
    pTable.age(neuronCounter) = expTable.Age(i);
    pTable.genotype(neuronCounter) = expTable.Genotype(i);
    pTable.GFPLabel(neuronCounter) = expTable.GFPLabel(i);
    pTable.eye(neuronCounter) = expTable.Eye(i);
    pTable.quad(neuronCounter) = expTable.Quadrant(i);
    pTable.SACside(neuronCounter) = expTable.SACside(i);
    pTable.stimSide(neuronCounter) = expTable.StimSide(i);
    pTable.cellNumb(neuronCounter) = expTable.CellNumb(i);
    pTable.distance(neuronCounter) = expTable.Distance(i);
    pTable.overlap(neuronCounter) = expTable.overlap(i);
    pTable.holdVol(neuronCounter) = expTable.holdVol(i);
    pTable.VhCorr{neuronCounter} = VhCorrList; % Store as cell array
    pTable.peakCurr{neuronCounter} = peakCurrList; % Store as cell array
    pTable.peakCurrCorr{neuronCounter} = peakCurrCorrList; % Store as cell array
    pTable.C2corrected{neuronCounter} = C2correctedList; % Store as cell array
    pTable.C1norm(neuronCounter, :) = C1norm; % Store normalized C1 data
    pTable.Rs{neuronCounter} = RsList;
    pTable.Rin{neuronCounter} = RinList;
    pTable.totalR{neuronCounter} = totalRList;

    % Increment the neuron counter
    neuronCounter = neuronCounter + 1;
end

% Save the final pTable
% save('pTable.mat', 'pTable');

