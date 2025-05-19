%% Build the neuron table FOR CURRENT CLAMP DATA 

%This code loads the data into iTable WITH PLOTS. See commented sections below 
% This code needs an excel spreadsheet with metadata in it.

% Every row in that excel spreadsheet is a single recording
% The buildTable code will run all the recording
%
% The end result of this code is that it saves 'iTable.mat'
% This is the table where each row is a recording
% Other codes will load this table and run analysis on it
%necessary functions are at the bottom of this script

%% Load the table
data_guide_name = '';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

 
% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'ExperimentDate','macPath','path','FileName','Sex','Age', 'Genotype', 'GFPLabel',...
    'Eye', 'Quadrant','RecordingType','Stim','StimFile','RMP','CellNumber','Drugs','barLength','barSpeed','numberOfRecs'},'string');

expTable = readtable(data_guide_name, opts');

[num_files, dummy] = size(expTable); %determines number of rows (iterations for loop will go through)
 

%% Intialize the neuron table (called iTable)

totalRecs = height(expTable);
neuronCounter = 1; %This value will keep track of which neuron I am on

dTable = table('Size', [totalRecs 16], 'VariableTypes', {'double','double','string','string','string',...
    'string','string','string','string','string','string','string','string','double','string','double'});

dTable.Properties.VariableNames = {'expDate', 'aquisitionNumber', 'sex','age', 'genotype',...
    'GFPLabel', 'eye', 'quad','recType','holdVol','stim','stimFile','stimTimes','cellNumb','drugs','holdCur'};

%% For loop for going through each recording to populate table

for i = 1:num_files
%this will run through all files
% for i = 1:2 % For testing code, num_files for all files
    
% changing the directory to go to the folder where the data is according to date 
    cd (expTable.path(i))
    % cd (expTable.macPath(i))

%creating variables to identify these properties for figure titles
    cellNumb = expTable.CellNumber(i);
    GFPlabel= expTable.GFPLabel(i);
    expDate = expTable.ExperimentDate(i);
    holdVol = expTable.holdVol(i);
    drugCond = expTable.Drugs(i);

    %for bars
    if expTable.Stim(i) == "bars"
    %load abf file
       [d,si]=abfload(expTable.FileName(i));

        %normalize baseline 
            x = mean(d(1:3000,:,:)); %averaging across the first .3 seconds 
            d_corr = d(:,:,:) - x;
            d =squeeze(d_corr);

            % Define filter parameters
            cutoff_frequency = 30; % Cut-off frequency in Hz (adjust as needed)
            sampling_rate = 10000; % Sample rate in Hz (adjust as needed)

            % Calculate the sampling interval
            dt = 1 / sampling_rate;

            % Apply the lowpass filter with IIR design
%             dFiltered = lowpass(d, cutoff_frequency, sampling_rate, 'ImpulseResponse', 'iir', 'Steepness', 0.5);
                dFiltered = lowpass(d, cutoff_frequency * 0.8, sampling_rate, 'ImpulseResponse', 'iir', 'Steepness', 0.9);

% %plots filtered raw data, comment to suppress raw data plots 
       hF = rawTest(d,si,cellNumb,expDate);

%% Loading stim files, directions, and times + orienting stim based on left and right eye 

       if expTable.AquisitionNumber(i) < 10 
         quickStim(append("00",num2str(expTable.AquisitionNumber(i))));
       end 
       if expTable.AquisitionNumber(i) > 9 
         quickStim(append("0",num2str(expTable.AquisitionNumber(i))));  
       end 
    

 %selects just the stimulus direction which is the first column 
       stimDir = stimInfo(1:24,1);
        Eye = expTable.Eye(i);

       [newInd] = eyeCalSOS(stimDir,Eye); 

       stimDirOld = stimDir;
       stimDir(:,1) = newInd; 
% % %spikes - separate ON and OFF component
% calculate spike count
    if expTable.Drugs(i) == "none"
        thresh = -7; %default -7
        manual = 0; 
        dON = d(15000:25000,1:24);
        dOFF = d(25000:45000,1:24);

        %finding the peaks of the subthreshold vm 
        dFiltON = dFiltered(15000:25000,1:24);
        dFiltOFF = dFiltered(30000:45000,1:24);

        pSubThreshON = max(dFiltON);
        pSubthreshOFF = max(dFiltOFF);

        %spike counts for ON 
        [nPts,nTrls] = size(dON);
        dVector = reshape(dON,numel(dON),1);
        dFilter = filterData(dVector,si);
        thresh = abs(thresh);
        spTmsRawON = getSpikeTimesCurrentClamp(dFilter,si,'thresh',thresh,'manual',manual);
        secTrial = nPts * si * 1e-6; %convert to seconds
        spTmsON = cell(nTrls,1);
        spCtsON = NaN(nTrls,1);

        for f = 1:nTrls
            spIndx = spTmsRawON > (f-1)*secTrial & spTmsRawON < f*secTrial;
            spikesON = spTmsRawON(spIndx) - (f-1)*secTrial;
            spTmsON{f} = spikesON;
            spCtsON(f) = numel(spikesON);
        end
        assignin('caller','spTmsON',spTmsON);
        assignin('caller','spCtsON',spCtsON);

        %spike counts for OFF
        [nPts,nTrls] = size(dOFF);
        dVector = reshape(dOFF,numel(dOFF),1);
        dFilter = filterData(dVector,si);
        thresh = abs(thresh);
        spTmsRawOFF = getSpikeTimesCurrentClamp(dFilter,si,'thresh',thresh,'manual',manual);
        secTrial = nPts * si * 1e-6; %convert to seconds
        spTmsOFF= cell(nTrls,1);
        spCtsOFF = NaN(nTrls,1);

        for f = 1:nTrls
            spIndx = spTmsRawOFF > (f-1)*secTrial & spTmsRawOFF < f*secTrial;
            spikesOFF = spTmsRawOFF(spIndx) - (f-1)*secTrial;
            spTmsOFF{f} = spikesOFF;
            spCtsOFF(f) = numel(spikesOFF);
        end
        assignin('caller','spTmsOFF',spTmsOFF);
        assignin('caller','spCtsOFF',spCtsOFF);

        % Sort data for ON and OFF based on direction & plot them
        pSortON = quickSort(spCtsON,stimDir);
        pSortOFF = quickSort(spCtsOFF,stimDir);
        %returns 8x3 double with each row is spike counts for the 8
        %different directions and each column is a trial 
        dirList = [0; 45; 90; 135; 180; 225; 270; 315];

        %do the same but for subthresh
        pFiltSortON = quickSort(pSubThreshON,stimDir);
        pFiltSortOFF = quickSort(pSubthreshOFF,stimDir);
        dirList = [0; 45; 90; 135; 180; 225; 270; 315];


% %% get DSI, vec length and pref dir for ON and OFF - uncomment in function to plot
        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSortON,stimDir,cellNumb,expDate,holdVol,drugCond,true);

        prefDirOn = prefDir; 
        DSIon =DSI; 
        prefSpikesON = prefSpikes; 
        nullSpikesON = nullSpikes; 
        %firing rate spikes/second based on the 15000 time window
        PrefFRon = prefSpikesON/15;
        NullFRon = nullSpikesON/15;

        if prefDirOn > 337.5
            prefDirOn = 0;
        end 

        [temp, prefInd] = min(abs(dirList - prefDir));
        absPrefDirOn = dirList(prefInd); 
        absNullDirOn = absPrefDirOn - 180;
        if absNullDirOn < 0 
            absNullDirOn = absNullDirOn + 360;
        end 

        %for OFF
        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSortOFF,stimDir,cellNumb,expDate,holdVol,drugCond,true);

        prefDirOff = prefDir; 
        DSIoff =DSI; 
        prefSpikesOFF = prefSpikes; 
        nullSpikesOFF = nullSpikes; 
        PrefFRoff = prefSpikesOFF/13; %firing rate spikes/second based on the 15000 time window
        NullFRoff = nullSpikesOFF/13;

        if prefDirOff > 337.5
            prefDirOff = 0;
        end 

        [temp, prefInd] = min(abs(dirList - prefDir));
        absPrefDirOff = dirList(prefInd); 
        absNullDirOff= absPrefDirOff - 180;

        if absNullDirOff < 0 
            absNullDirOff = absNullDirOff + 360;
        end 

        [avgTrace,prefDirTrace,nullDirTrace] = plotDirTraces(d,stimDir,1,10000,absNullDirOn,absPrefDirOn);

        SpikeprefDirTrace = mean(prefDirTrace,2);
        SpikenullDirTrace = mean(nullDirTrace,2);

        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pFiltSortON,stimDir,cellNumb,expDate,holdVol,drugCond,true);
        SubprefDirOn = prefDir; 
        SubDSIon = abs(DSI);
        SubVecSumOn = vecLength;

        if SubprefDirOn > 337.5
            SubprefDirOn = 0;
        end 

        [temp, prefInd] = min(abs(dirList - prefDir));
        SubabsPrefDirOn = dirList(prefInd); 
        SubabsNullDirOn = SubabsPrefDirOn - 180;
        if SubabsNullDirOn < 0 
            SubabsNullDirOn = SubabsNullDirOn + 360;
        end 
% 
        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pFiltSortOFF,stimDir,cellNumb,expDate,holdVol,drugCond,true);
        SubprefDirOff = prefDir; 
        SubDSIoff = abs(DSI);
        SubVecSumOff = vecLength;

        if SubprefDirOn > 337.5
            SubprefDirOn = 0;
        end 

        [temp, prefInd] = min(abs(dirList - prefDir));
        SubabsPrefDirOff = dirList(prefInd); 
        SubabsNullDirOff = SubabsPrefDirOff - 180;
        if SubabsNullDirOff < 0 
            SubabsNullDirOff = SubabsNullDirOff + 360;
        end 

        [avgTrace,prefDirTrace,nullDirTrace] = plotDirTraces(dFiltered,stimDir,1,10000,absNullDirOn,absPrefDirOn);

        SubprefDirTrace = mean(prefDirTrace,2);
        SubnullDirTrace = mean(nullDirTrace,2);


        %norm VS for sub
            %for ON
            pMeanON = mean(pFiltSortON,2); %average response for each direction 
            pMaxON = max(pMeanON); %finds max response
            pNormON = pMeanON/pMaxON; %normalizes all resp to max resp
            %for OFF
            pMeanOFF = mean(pFiltSortOFF,2); 
            pMaxOFF = max(pMeanOFF); 
            pNormOFF = pMeanOFF/pMaxOFF;

            nReps = size(pFiltSortOFF,2); %defining number of repitions 
            uDirs = unique(stimDir);% list of directions 
            uDirs = deg2rad(uDirs); %converts dirs to radians 

            %normVS for ON 
            [x,y] = pol2cart(uDirs, pNormON); %polar to cartesian 
            prefDir = atan2d(sum(y),sum(x)); %vector sum 
            if prefDir < 0
                  prefDir = prefDir + 360;
            end
            normVSon = sqrt(sum(x)^2 + sum(y)^2);

            %norm for OFF
            [x,y] = pol2cart(uDirs, pNormOFF); %polar to cartesian 
            prefDir = atan2d(sum(y),sum(x)); %vector sum 
            if prefDir < 0
                  prefDir = prefDir + 360;
            end
            normVSoff = sqrt(sum(x)^2 + sum(y)^2);


    %% get average peak response for pref and null dir
    x = stimDir; 
    y = pSubThreshON'; 
    z = pSubthreshOFF';

    %for pref ON
    j = find(x==absPrefDirOn);
    k = y(j); 
    pSubPrefON = mean(k); 
    %for null ON
    l = find(x==absNullDirOn);
    m =y(l);
    pSubNullON = mean(m);
    %for pref OFF
    n = find(x==absPrefDirOff);
    o = z(n);
    pSubPrefOFF = mean(o);
    %for null OFF
    p = find(x==absNullDirOff);
    q = z(p);
    pSubNullOFF = mean(q); 

%% %% Maximum firing rate

%FOR ON 
maxFiringRatesON = cell(24, 1);
maxFiringRatesOFF = cell(24, 1);

for x = 1:24
    spCts = spCtsON(x,:); 
    spCts = ones(x,spCts)';
    spTms = spTmsON{x,:}';
    spikenTimes = [spCts spTms];
    % Define the parameters
    windowSize = 0.1; % Size of the sliding window in seconds
    totalTime = 1.5; % Total duration of the spike recording
    numWindows = floor(totalTime / windowSize); % Calculate the number of windows

    % Initialize variables to store firing rates and window start times
    firingRates = zeros(numWindows, 1);
    windowStartTimes = (0:numWindows-1) * windowSize;

% Slide the window across the spike time data
for j = 1:numWindows
    startTime = (j - 1) * windowSize;
    endTime = startTime + windowSize;

    % Count the number of spikes within the current window
    spikesInWindow = sum(spikenTimes(:, 2) >= startTime & spikenTimes(:, 2) < endTime);

    % Calculate the firing rate for the current window
    firingRates(j) = spikesInWindow / windowSize;
end

% Find the window with the maximum firing rate
[maxFiringRate, maxIndex] = max(firingRates);
maxStartTime = windowStartTimes(maxIndex);

maxFiringRatesON{x} = max(firingRates);
prefDirIndices = find(stimDir == absPrefDirOn);
nullDirIndices = find(stimDir == absNullDirOn);
    for y = 1:numel(prefDirIndices)
        index = prefDirIndices(y);
        firingRate = maxFiringRatesON{index};
        firingRatesPrefDir{y} = firingRate;
    end
    for y = 1:numel(nullDirIndices)
        index = nullDirIndices(y);
        firingRate = maxFiringRatesON{index};
        firingRatesNullDir{y} = firingRate;
    end
    avgFRprefOn = mean(cell2mat(firingRatesPrefDir));
    avgFRnullOn = mean(cell2mat(firingRatesNullDir));
end

% for OFF
for x = 1:24
    spCts = spCtsOFF(x,:); 
    spCts = ones(x,spCts)';
    spTms = spTmsOFF{x,:}';
    spikenTimes = [spCts spTms];
% Define the parameters
windowSize = 0.1; % Size of the sliding window in seconds
totalTime = 1.5; % Total duration of the spike recording
numWindows = floor(totalTime / windowSize); % Calculate the number of windows

% Initialize variables to store firing rates and window start times
firingRates = zeros(numWindows, 1);
windowStartTimes = (0:numWindows-1) * windowSize;

% Slide the window across the spike time data
for j = 1:numWindows
    startTime = (j - 1) * windowSize;
    endTime = startTime + windowSize;

    % Count the number of spikes within the current window
    spikesInWindow = sum(spikenTimes(:, 2) >= startTime & spikenTimes(:, 2) < endTime);

    % Calculate the firing rate for the current window
    firingRates(j) = spikesInWindow / windowSize;
end

% Find the window with the maximum firing rate
[maxFiringRate, maxIndex] = max(firingRates);
maxStartTime = windowStartTimes(maxIndex);

maxFiringRatesOFF{x} = max(firingRates);
prefDirIndices = find(stimDir == absPrefDirOff);
nullDirIndices = find(stimDir == absNullDirOff);
    for y = 1:numel(prefDirIndices)
        index = prefDirIndices(y);
        firingRate = maxFiringRatesOFF{index};
        firingRatesPrefDir{y} = firingRate;
    end
    for y = 1:numel(nullDirIndices)
        index = nullDirIndices(y);
        firingRate = maxFiringRatesOFF{index};
        firingRatesNullDir{y} = firingRate;
    end
    avgFRprefOff = mean(cell2mat(firingRatesPrefDir));
    avgFRnullOff = mean(cell2mat(firingRatesNullDir));
end 
    end 

    if expTable.Drugs(i) == "TTX"

        dON = d(15000:30000,1:24);
        dOFF = d(32000:45000,1:24);
        %finding the peaks of the subthreshold vm 
        pSubThreshON = max(dON);
        pSubthreshOFF = max(dOFF);
        %do the same but for subthresh
        pSubSortON = quickSort(pSubThreshON,stimDir);
        pSubSortOFF = quickSort(pSubthreshOFF,stimDir);
        dirList = [0; 45; 90; 135; 180; 225; 270; 315];

        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSubSortON, stimDir,cellNumb,expDate,holdVol,drugCond,true);
        SubprefDirOn = prefDir; 
        SubDSIon = abs(DSI);
        SubVecSumOn = vecLength;

        if SubprefDirOn > 337.5
            SubprefDirOn = 0;
        end 

        [temp, prefInd] = min(abs(dirList - prefDir));
        SubabsPrefDirOn = dirList(prefInd); 
        SubabsNullDirOn = SubabsPrefDirOn - 180;
        if SubabsNullDirOn < 0 
            SubabsNullDirOn = SubabsNullDirOn + 360;
        end 
% 
        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSubSortOFF,stimDir,cellNumb,expDate,holdVol,drugCond,true);
        SubprefDirOff = prefDir; 
        SubDSIoff = abs(DSI);
        SubVecSumOff = vecLength;

        if SubprefDirOn > 337.5
            SubprefDirOn = 0;
        end 

        [temp, prefInd] = min(abs(dirList - prefDir));
        SubabsPrefDirOff = dirList(prefInd); 
        SubabsNullDirOff = SubabsPrefDirOff - 180;
        if SubabsNullDirOff < 0 
            SubabsNullDirOff = SubabsNullDirOff + 360;
        end 

        [avgTrace,prefDirTrace,nullDirTrace] = plotDirTraces(d,stimDir,1,10000,SubabsNullDirOn,SubabsPrefDirOn);

        SubprefDirTrace = mean(prefDirTrace,2);
        SubnullDirTrace = mean(nullDirTrace,2);
    %norm VS for sub
            %for ON
            pMeanON = mean(pSubSortON,2); %average response for each direction 
            pMaxON = max(pMeanON); %finds max response
            pNormON = pMeanON/pMaxON; %normalizes all resp to max resp
            %for OFF
            pMeanOFF = mean(pSubSortOFF,2); 
            pMaxOFF = max(pMeanOFF); 
            pNormOFF = pMeanOFF/pMaxOFF;

            nReps = size(pSubSortOFF,2); %defining number of repitions 
            uDirs = unique(stimDir);% list of directions 
            uDirs = deg2rad(uDirs); %converts dirs to radians 

            %normVS for ON 
            [x,y] = pol2cart(uDirs, pNormON); %polar to cartesian 
            prefDir = atan2d(sum(y),sum(x)); %vector sum 
            if prefDir < 0
                  prefDir = prefDir + 360;
            end
            normVSon = sqrt(sum(x)^2 + sum(y)^2);

            %norm for OFF
            [x,y] = pol2cart(uDirs, pNormOFF); %polar to cartesian 
            prefDir = atan2d(sum(y),sum(x)); %vector sum 
            if prefDir < 0
                  prefDir = prefDir + 360;
            end
            normVSoff = sqrt(sum(x)^2 + sum(y)^2); 


    %% get average peak response for pref and null dir
    x = stimDir; 
    y = pSubThreshON'; 
    z = pSubthreshOFF';

    %for pref ON
    j = find(x==SubabsPrefDirOn);
    k = y(j); 
    pSubPrefON = mean(k); 
    %for null ON
    l = find(x==SubabsNullDirOn);
    m =y(l);
    pSubNullON = mean(m);
    %for pref OFF
    n = find(x==SubabsPrefDirOff);
    o = z(n);
    pSubPrefOFF = mean(o);
    %for null OFF
    p = find(x==SubabsNullDirOff);
    q = z(p);
    pSubNullOFF = mean(q); 


    end 



%% Load neuron infor into the dTable

    dTable.expDate(i) = expTable.ExperimentDate(i);
    dTable.aquisitionNumber(i) = expTable.AquisitionNumber(i);
    dTable.sex(i) = expTable.Sex(i);
    dTable.age(i) = expTable.Age(i);
    dTable.genotype(i) = expTable.Genotype(i);
    dTable.GFPLabel(i) = expTable.GFPLabel(i);
    dTable.eye(i) = expTable.Eye(i);
    dTable.quad(i) = expTable.Quadrant(i);
    dTable.recType(i) = expTable.RecordingType(i);
    dTable.holdVol(i) = expTable.holdVol(i);
    dTable.stim(i) = expTable.Stim(i);
    dTable.stimFile(i) = expTable.StimFile(i);
    dTable.cellNumb(i) = expTable.CellNumber(i);
    dTable.drugs(i) = expTable.Drugs(i);
    dTable.quad(i) = expTable.Quadrant(i);
    dTable.SubprefDirTrace{i} = SubprefDirTrace;
    dTable.SubnullDirTrace{i} = SubnullDirTrace;
    dTable.avgFRprefOn(i) = avgFRprefOn;
    dTable.avgFRnullOn(i) = avgFRnullOn;
    dTable.avgFRprefOff(i) = avgFRprefOff;
    dTable.avgFRnullOff(i) = avgFRnullOff;
    dTable.DSIon(i) = DSIon;
    dTable.DSIoff(i) = DSIoff;
    dTable.prefDirOn(i) = prefDirOn;
    dTable.prefDirOff(i) = prefDirOff;
    dTable.absPrefDirOn(i) = absPrefDirOn;
    dTable.absPrefDirOff(i) = absPrefDirOff;
    dTable.absNullDirOn(i) = absNullDirOn;
    dTable.absNullDirOff(i) = absNullDirOff;
    dTable.prefSpikesON(i) = prefSpikesON;
    dTable.prefSpikesOFF(i) = prefSpikesOFF;
    dTable.nullSpikesON(i)= nullSpikesON;
    dTable.nullSpikesOFF(i) = nullSpikesOFF;
    dTable.pSubPrefON(i) = pSubPrefON;
    dTable.pSubPrefOFF(i) = pSubPrefOFF;
    dTable.pSubNullON(i) = pSubNullON;
    dTable.pSubNullOFF(i) = pSubNullOFF;
    dTable.SubDSIon(i) = SubDSIon;
    dTable.SubDSIoff(i) = SubDSIoff;
    dTable.SubnormVSon(i) = normVSon; 
    dTable.SubnormVSoff(i) = normVSoff; 
    dTable.SubVecSumOn(i) = SubVecSumOn;
    dTable.SubVecSumOff(i) = SubVecSumOff;
    dTable.SubprefDirOn(i) = SubprefDirOn;
    dTable.SubprefDirOff(i) = SubprefDirOff;
    dTable.SubabsPrefDirOn(i) = SubabsPrefDirOn;
    dTable.SubabsPrefDirOff(i) = SubabsPrefDirOff;
    dTable.SubabsNullDirOn(i) = SubabsNullDirOn;
    dTable.SubabsNullDirOff(i) = SubabsNullDirOff;
    dTable.FRprefON(i) = PrefFRon;
    dTable.FRnullON(i) = NullFRon;
    dTable.FRprefOFF(i) = PrefFRoff;
    dTable.FRnullOFF(i) = NullFRoff;
    dTable.spikeCountOn(i,:) = spCtsON';
    dTable.spikeCountOFF(i,:) = spCtsOFF';

 end 
  end 
%     end 
% dTable = dTable(strcmp(dTable.drugs,"TTX"),:);
% dTable = dTable(strcmp(dTable.drugs,"none"),:);
dTable = dTable(strcmp(dTable.stim,"bars"),:);
%     end 
    
%% Functions 

function hF = rawTest(d,si,cellNumb,expDate)
%plots unprocessed data if called, outputs figure handle.
%nPts = number of points (time), trials = number of directions 

[nPts,trials] = size(d);
dt = si*1e-6;
t = 0:dt:(nPts - 1)*dt;
L = floor(sqrt(trials));
W = ceil(trials / L);
minD = min(d(:));
maxD = max(d(:));
hF = figure;

titleStr = sprintf('Neuron %g Date %g',cellNumb,expDate);

for i = 1:trials
    subplot(L,W,i)
    plot(t,d(:,i),'r') 
    axis tight
    ylim([minD 100])
    
end
title(titleStr);
end

function quickStim(stimNum,stimDate)
% Mathew's function to load stimXXX.txt files for a given experiment.

if nargin < 2 || isempty(stimDate) 
    %if no input argument, assume current directory.
    [~,stimDate] = fileparts(pwd);
    fprintf('Loading from current folder; %s \n',stimDate);
else
    stimDate = char(stimDate); %ensure indexable characters, not a string
end

dsgcDir = '/Users/karina bistrong/Dropbox/Mac/Desktop/DATA/Ephys';
searchDirName = sprintf('%s*',stimDate); %find directories that match input date

oldDir = cd;
newDir = dir(searchDirName);
cd(newDir.name);

fn = sprintf('stim%s.txt',stimNum);
stimInfo = load('-ASCII',fn);
assignin('caller','stimInfo',stimInfo);
% order of stimInfo depends on stim function used, but usually:
% 1st column: bar directions
% 2nd column: bar speeds
% 3rd column: bar length
% 4th column: bar width
% 5th column: radius to traverse
% 6th column: bar brightness

cd(oldDir);

end

function ySort = quickSort(y,x)
%Sorts a vector of elements y (e.g. spike counts, or spike times) based on
%the unique elements of an input vector "x" (e.g. stim conditions)

uX = unique(x);
nuX = numel(uX); %number of elements 
[nDim1,nDim2] = size(y);

if nDim2 == 1 %if dim 2 is singleton, properly align dimensions
    nPts = nDim2;
    nTrials = nDim1;
else %otherwise assume dim1 is nPoints and dim2 nTrials
    nPts = nDim1;
    nTrials = nDim2;
end

nReps = nTrials / nuX;
sortFlag = 0;
if iscell(y)
    ySort = cell(nuX,nReps);
elseif nPts > 1
    ySort = zeros(nPts,nReps,nuX);
    sortFlag = 1;
else
    ySort = zeros(nuX,nReps);
end

if sortFlag
    
    for i = 1:nuX
        indx = ( uX(i) == x ); %find each index corresponding to a given stim
        ySort(:,:,i) = y(:,indx);
    end
    
else
    
    for i = 1:nuX
        indx = ( uX(i) == x ); %find each index corresponding to a given stim
        ySort(i,:) = y(indx);
    end
    
end 
end
% 
function [avgTrace, prefDirTrace, nullDirTrace] = plotDirTraces(stimDF,stimDir,ROI,Fs,absNullDirOn,absPrefDirOn)

if nargin < 4 || isempty(Fs) %if no Fs assume Ca2+ imaging data
    Fs = 1.48;
    LW = 1.5;
    scaleBar = [0 .2];
elseif Fs == 1.48 || Fs == 2.96 %Ca2+ imaging data
    LW = 1.5;
    scaleBar = [0 .2];
else %otherwise assume ephys
    LW = .5;
    if sum(stimDF(:)) > 0 %try to determine if primarily positive or negative signal
        scaleBar = [0 20];
    else
        scaleBar = [-20 0];
    end
end

dFSort = quickSort(stimDF(:,:,ROI),stimDir); %sorts the stimDF data by the stimulus direction 

[nFrames,nReps,nDirs] = size(dFSort);
nTrials = nReps*nDirs;
uDirs = unique(stimDir);
plotDF = reshape(dFSort,nFrames,nTrials); %reshapes the data into a 2D matrix size of frames x trials where trials is the repsc

tFrames = 1:nFrames;
tSec = (tFrames - 1)/Fs;

mindF = min(dFSort(:));
maxdF = max(dFSort(:));

xStart = .05;
xEnd = .95;
xEach = (xEnd - xStart)/nDirs;

yStart = .05;
yEnd = .95;
yEach = (yEnd - yStart)/nReps;

yLabelWidth = 1 - yEnd;

hF = figure;
%co = hsv(nTrials);
for i = 1:nTrials
    xInt = floor((i-1)/nReps);%mod(i-1,nDirs);
    yInt = mod(i-1,nReps) + 1;%floor((i-1)/nDirs) + 1;
    %cInt = mod(i,2) + 1;
    pos = [(xStart + xInt*xEach) (yEnd - yInt*yEach) xEach yEach];
    hA = axes('Units','Normalized','Position',pos,'XTick',[],'YTick',[],...
        'YLim',[mindF 100],'XLim',[tSec(1) tSec(end)],'NextPlot','replacechildren');
    plot(tSec,plotDF(:,i),'k','LineWidth',LW)
end

hF.Children(end).YTick = scaleBar;
%hF.Children(end).YLabel.String = 'deltaF/F';

for j = 1:nDirs
    txtPos = [(xStart + (j-1)*xEach) yEnd xEach yLabelWidth];
    txtStr = num2str(uDirs(j));
    txt = uicontrol(hF,'Style','text','Units','Normalized','Position',txtPos,...
        'String',txtStr,'FontSize',12);
end
%hA.YAxisLocation = 'right';
%hA.YTick = [0 (maxdF - mindF)/2 ];



hF=figure; 
for j = 1:nDirs
    % plot average trace
    subplot(nReps+1, nDirs, nReps*nDirs+j)
    avgTrace = mean(dFSort(:, 1:nReps, j), 2);
    plot(tSec, avgTrace, 'k', 'LineWidth', LW)
%     xlabel('Time (s)')
%     ylabel('\DeltaF/F')
    ylim([mindF 100])
    xlim([tSec(1) tSec(end)])
    title([ num2str(uDirs(j))])
    txtPos = [(xStart + (j-1)*xEach) yEnd xEach yLabelWidth];
    txtStr = num2str(uDirs(j));
    txt = uicontrol(hF,'Style','text','Units','Normalized','Position',txtPos,...
        'String',txtStr,'FontSize',12);
end 
% Extract the subplot for the preferred direction "0"
prefDirIndex = find(uDirs == absPrefDirOn); % Find the index of the preferred direction "180"
prefDirTrace = dFSort(:, 1:nReps, prefDirIndex); % Extract the data for the preferred direction
% Extract the subplot for the null direction "180"
nullDirIndex = find(uDirs == absNullDirOn); % Find the index of the null direction "0"
nullDirTrace = dFSort(:, 1:nReps, nullDirIndex); % Extract the data for the null direction

end

% 
function [hF,prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSort,stimDir,cellNumb,expDate,holdVol,drugCond,showLess,makeFig)

 hF = [];
if nargin<8|| isempty(makeFig)
    makeFig = 1;
end
if nargin<7 || isempty(showLess)
    showLess = 0;
end

ctMean = mean(pSort,2); %takes the mean of each row in the sorted ON or OFF responses = average response for each direction 
ctSortPlot = [pSort; pSort(1,:)]; %produces a matrix of the sorted responses ie is just pSortON or pSort OFF
ctMeanPlot = [ctMean; ctMean(1,:)]; %matrix of the means 

nReps = size(pSort,2); %defining number of repitions of directions based on number of columns in pSort 

uDirs = unique(stimDir);% list of directions %potentially problematic way of doing things;
% if stimDirs is shifted in a non-uniform way (e.g. resets to 0 after passing 360)
% then ctSort will now be incorrectly indexed
uDirs = deg2rad(uDirs); %converts list of directions to radians 

[x,y] = pol2cart(uDirs, ctMean / sum(ctMean)); %polar coordinates to cartesian coordinates
prefDir = atan2d(sum(y),sum(x)); %vector sum 
if prefDir < 0
    prefDir = prefDir + 360;
end
vecLength = sqrt(sum(x)^2 + sum(y)^2);

nDirs = length(uDirs);  %number of directions = 8 
incDirs = 360 / nDirs; %direction increments 
prefIndx = round(prefDir / incDirs) + 1;
if prefIndx > nDirs
    prefIndx = 1;
end
nullIndx = prefIndx - (nDirs / 2);
if nullIndx < 1
    nullIndx = nullIndx + nDirs;
end
DSI = (ctMean(prefIndx) - ctMean(nullIndx) ) / (ctMean(prefIndx) + ctMean(nullIndx));
preSet_rlim = [0 40];
prefSpikes = ctMean(prefIndx);
nullSpikes = ctMean(nullIndx);
titleStr = sprintf('Pref Dir %3.1f DSI %4.2f Vec Length %4.2f',prefDir,DSI,vecLength);
% acqInfo = sprintf('Neuron %g Date %g Hold %g Drug %s',cellNumb,expDate,holdVol,drugCond);
uDirs = [uDirs; uDirs(1)];
uDirsPlot = repmat(uDirs,1,nReps);


% Plot figure
if ~showLess
    if makeFig
        hF = figure;
    else
        hF = [];
    end

    %Set radial bounds
    pAx = polaraxes();
    
    tiledlayout(2,1);
   

    plotChild = polarplot(uDirsPlot,ctSortPlot);
    set(plotChild(1:nReps),'LineWidth',1);
    
    hold on
 
    plotChild = polarplot(uDirs,ctMeanPlot,'k');
    set(plotChild(1),'LineWidth',2)
    
    
%     plotChild = compass(sum(x)*max(ctMean),sum(y)*max(ctMean),'k'); %think about this
    plotChild = polarplot(deg2rad([prefDir, prefDir]), [0 max(ctMean)*vecLength],'k');
    set(plotChild(1),'LineWidth',1.5)
    
    % title(titleStr,acqInfo);
    rlim(preSet_rlim);
%     if rFlag
%         pAx.RLim = rBounds;
%     end
    
% else
    hF = [];
end
end
function [spike_times,excluded_trials,fig_num] = getSpikeTimesCurrentClamp(data,samp_int,varargin)
% finds spiketimes (local minima) in clampex data loaded with abfload
%
% %INPUTS:
%     data                                      Nx1 vector
%     samp_int                            scalar sampling interval in microseconds
%     thresh                                 threshold in std
%     manual                                 0 or 1, manually choose spike threshold
%     show_plot                          0 or 1, optional plot flag
%     ax_offset                          time offset in samples for plot
%     y_axis_range                  'auto' or [min max]
%     refractory_period     scalar, time in ms to exclude spikes
%     tlt                                          titlr of the plot
%     dirs                                        if a whole repetition that contains all different directions is given, 
%                                                       the plot will be split by lines between the trials of the different dirs
%     dur                                          the duration of each trial (in sec)
%
% %OUTPUTS
%     SpikeTimes                        1xN vector of spiketimes
%     excluded_trials           the nmbers of all trials that should be excluded from the analysis
%     fig_num                                the figure number, so we can close it later

%default values
thresh = 7;
manual = 0;
show_plot=0;
ax_offset = 0;
y_axis_range = 'auto';
refractory_period = 0.001;
tlt = '';
dirs = [];
dur=0;

%assign input variables
pvpmod(varargin);

default_cutoff=mean(data)+thresh*std(data);
spike_times = [];
excluded_trials = [];
fig_num = [];
lower_cutoff = 0;

if manual
    [~,~,fig_num] = makeTracePlot(data,samp_int,'ax_offset',ax_offset,'y_axis_range',y_axis_range);
    [~,ys] = selectThreshold(default_cutoff);
    if ~isempty(ys)
        cutoff = max(ys);
        lower_cutoff = min(ys);
    else %no spikes deined
        cutoff = min(data)-1;
    end
    close;
else
    cutoff=default_cutoff;    
end

%find spikes according to local minima
ind = 1;
suspect_points = find(data>=cutoff);    
for d = 2:(length(suspect_points)-1)
    if data(suspect_points(d)) >= data(suspect_points(d)-1) && data(suspect_points(d)) >= data(suspect_points(d)+1)
        spike_times(ind) = (suspect_points(d))*samp_int*1e-6; %convert to seconds
        spike_vals(ind) = data(suspect_points(d));
        ind = ind + 1;
    end
end
%remove refractory period violations
bad_inds = diff(spike_times) <= refractory_period;
spike_times(bad_inds) = [];
spike_vals(bad_inds) = [];

% remove noises with a high amplitude
if lower_cutoff
    rem = find(spike_vals<lower_cutoff);
    spike_times(rem) = [];
    spike_vals(rem) = [];
end

%convert sec
ax_offset_sec = ax_offset*samp_int*1e-6; 
spike_times = spike_times + ax_offset_sec;

%optional plot of data with selected threshold
if show_plot
    [~,ax,fig_num] = makeTracePlot(data,samp_int,'ax_offset',ax_offset);
    %show threshold cutoff for spikes
    hold on
    title(sprintf('%s; %d spikes',tlt,length(spike_times)));
    threshline = ones(length(data),1).*cutoff;
    plot(ax,threshline,'r--','LineWidth',2);
    if spike_times
        plot(spike_times,spike_vals,'r.','MarkerSize',10);
    end
    hold off
    ylim(y_axis_range);
    if ~isempty(dirs)
        trial_list = cell(length(dirs),1);
        for d=1:length(dirs)
            trial_list{d} = num2str(d);
            line(dur*(d-1)*[1 1],[min(data) 1.2*max(data)],'color','r');
% % %             % Temporal change: %%%%%%%
% % %             f=sprintf('text(%d, %d, ''%d'');',dur*(d-1),1.1*max(data),dirs(d));            
            % Do not write the directions, so no bias in the decision of trials to exclude            
            % Instead, just give each trial a number for recognition
            f=sprintf('text(%d, %d, ''%d'');',dur*(d-1),1.1*max(data),d);
            eval(f);
        end
        if manual
            [excluded_trials,~] = listdlg('ListString',trial_list,'CancelString','Include All','PromptString','Select trials to exclude');
        end
    end
end

return
end 