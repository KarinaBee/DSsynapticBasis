%% Build the neuron table FOR CELL ATTACHED DATA 

%This code just loads the data into cTable WITH PLOTS. See commented sections below 

% This code needs an excel spreadsheet with metadata in it.

% Every row in that excel spreadsheet is a single recording
% The buildTable code will run all the recording
%
% The end result of this code is that it saves 'cTable.mat'
% This is the table where each row is a recording
% Other codes will load this table and run analysis on it
%necessary functions are at the bottom of the script 

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
    barSpeed = expTable.barSpeed(i);
        %load abf file
       [d,si]=abfload(expTable.FileName(i));

       %normalize baseline 
            x = mean(d(1:3000,:,:)); %averaging across the first .3 seconds 
            d_corr = d(:,:,:) - x;
            d =squeeze(d_corr);
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
       newInd = zeros(24,1); %new empty list for new stim dirs 
       temp = unique(stimDir(:,1));   

       %Right eye calibration
       right0 =180;
       right45 =135;
       right90 =90;
       right135 =45; 
       right180 =0;
       right225 =315;
       right270 =270;
       right315 =225; 

       if expTable.Eye(i) == "R"
        replaceDir = [right0; right45; right90; right135; right180; right225; right270; right315];
       end 

       %Left/ventral eye calibration 
       left0 = 0; 
       left45 = 45;
       left90 = 90; 
       left135 = 135;
       left180 = 180;
       left225 = 225;
       left270 = 270;
       left315 = 315; 

       if expTable.Eye(i) == "L"
        replaceDir = [left0; left45; left90; left135; left180; left225; left270; left315];
       end 

       for ii = 1:numel(temp)
        indices = stimDir(:,1) == temp(ii);
        newInd(indices) = replaceDir(ii);
       end

       stimDirOld = stimDir;
       stimDir = newInd;

        % Parameters
        thresh = -7; % Default threshold value for spike detection
        manual = 0; % Manual flag

        % Data segmentation for ON and OFF phases
        if expTable.barSpeed(i) == "500"
        dON = d(15000:30000, 1:24);
        dOFF = d(30000:45000, 1:24);
        end 

     
        if expTable.barSpeed(i) == "250"
        dON = d(35000:50000, 1:24);
        dOFF = d(55000:70000, 1:24);
        end 

        % Plot filtered data for ON and OFF phases
        figure;
        subplot(2, 1, 1);
        plot(dON);
        title('ON Phase');
        subplot(2, 1, 2);
        plot(dOFF);
        title('OFF Phase');

        % Process ON phase
        [spTmsON, spCtsON] = process_phase(dON, thresh, si, manual);
        assignin('caller', 'spTmsON', spTmsON);
        assignin('caller', 'spCtsON', spCtsON);

        % Process OFF phase
        [spTmsOFF, spCtsOFF] = process_phase(dOFF, thresh, si, manual);
        assignin('caller', 'spTmsOFF', spTmsOFF);
        assignin('caller', 'spCtsOFF', spCtsOFF);

      % Process entire trace
        [spTmsALL, spCtsALL] = process_phase(d, thresh, si, manual);
        assignin('caller', 'spTmsALL', spTmsALL);
        assignin('caller', 'spCtsALL', spCtsALL);

        pSortAll = quickSort(spCtsALL,stimDir);        

        pSortON = quickSort(spCtsON,stimDir);
        pSortOFF = quickSort(spCtsOFF,stimDir);
        %returns 8x3 double with each row is spike counts for the 8
        %different directions and each column is a trial 
        dirList = [0; 45; 90; 135; 180; 225; 270; 315];

        %get DSI, vec length and pref dir for ON and OFF - uncomment in function to plot
        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSortON,stimDir,cellNumb,expDate,holdVol,drugCond);

        prefDirOn = prefDir; 
        DSIon =DSI; 
        prefSpikesON = prefSpikes; 
        nullSpikesON = nullSpikes; 

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
        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSortOFF,stimDir,cellNumb,expDate,holdVol,drugCond);

        prefDirOff = prefDir; 
        DSIoff =DSI; 
        prefSpikesOFF = prefSpikes; 
        nullSpikesOFF = nullSpikes; 

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

        PDdiff = abs(absPrefDirOn - absPrefDirOff);

        %normVS for SpikeS ON 
          %for ON
            pMeanON = mean(pSortON,2); %average response for each direction 
            pMaxON = max(pMeanON); %finds max response
            pNormON = pMeanON/pMaxON; %normalizes all resp to max resp
            %for OFF
            pMeanOFF = mean(pSortOFF,2); 
            pMaxOFF = max(pMeanOFF); 
            pNormOFF = pMeanOFF/pMaxOFF;

            nReps = size(pSortOFF,2); %defining number of repitions 
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


%% %% Maximum firing rate

% Initialize variables
%for on and off
maxFRon = zeros(24,1);
maxFRoff = zeros(24,1);
firingRatesOnWindow = zeros(24,15);
firingRatesOffWindow = zeros(24,15);

%for whole trace
firingRatesAllWindows = zeros(24, 45);
allMaxFR = zeros(24,1); 

% Process each trial for whole trace 

for x = 1:24
    spCts = spCtsALL(x, :);
    spCts = ones(1, spCts)';
    spTms = spTmsALL{x, :}';
    spikeTimes = [spCts spTms];

    % Parameters for sliding window
     windowSize = 0.1; % 100ms sliding window
     totalTime = 4.5; % total duration of spike phase
     numWindows = floor(totalTime / windowSize);

    % Variables to store firing rates and window start times
    firingRates = zeros(numWindows, 1);
    windowStartTimes = (0:numWindows-1) * windowSize;

    % Sliding the window
    for j = 1:numWindows
        startTime = (j - 1) * windowSize;
        endTime = startTime + windowSize;

        % Count number of spikes in current window
        spikesInWindow = sum(spikeTimes(:,2) >= startTime & spikeTimes(:,2) < endTime);

        % Firing rate for current window
        firingRates(j) = spikesInWindow / windowSize;
    end

    % Store firing rates for each window for the current trial
    firingRatesAllWindows(x, :) = firingRates;

    % Find and store max firing rate for current trial
    allMaxFR(x) = max(firingRates);
end


% Initialize storage for preferred and null direction firing rates
allFiringRatesPD = [];
allFiringRatesND = [];

% Finding and calculating average FR for preferred direction (PD) ON
prefDirIndices = find(stimDir == absPrefDirOn);
for y = 1:numel(prefDirIndices)
    index = prefDirIndices(y);
    allFiringRatesPD = [allFiringRatesPD; firingRatesAllWindows(index, :)];
end


% Finding and calculating average FR for null direction (ND) ON
nullDirIndices = find(stimDir == absNullDirOn);
for y = 1:numel(nullDirIndices)
    index = nullDirIndices(y);
    allFiringRatesND = [allFiringRatesND; firingRatesAllWindows(index, :)];
end

% Calculate the average firing rates for preferred and null directions
avgPDfiringRate = mean(allFiringRatesPD, 1);
avgNDfiringRate = mean(allFiringRatesND, 1);

% Define the window bins
windowBins = (0:0.1:4.4); % 45 bins of size 0.1 seconds

%For ON

for x = 1:24
    spCts = spCtsON(x,:); 
    spCts = ones(1, spCts)'; 
    spTms = spTmsON{x,:}';
    spikenTimes = [spCts spTms];

    % Define parameters for sliding window 
    windowSize = 0.1; % 100ms sliding window
    totalTime = 1.5; % total duration of spike phase 
    numWindows = floor(totalTime / windowSize); 

    % Variables to store firing rates and window start times
    firingRates = zeros(numWindows, 1); 

    % Sliding the window 
    for j = 1:numWindows
        startTime = (j - 1) * windowSize; 
        endTime = startTime + windowSize; 

        % Count number of spikes in current window
        spikesInWindow = sum(spikenTimes(:,2) >= startTime & spikenTimes(:,2) < endTime);
        firingRates(j) = spikesInWindow / windowSize; % Firing rate for current window
    end 

    % Store firing rates for each window for the current trial
    firingRatesOnWindow(x, :) = firingRates; 
end

%  disp(maxFRon); 

    %Finding + caclulating average max FR for PD off 
    prefDirIndices = find(stimDir == absPrefDirOn);
    for y = 1:numel(prefDirIndices)
        index = prefDirIndices(y);
        maxFR = maxFRon(index);
        firingRatesPrefDir{y} = maxFR;
    end
%     disp(firingRatesPrefDir);
    avgFRprefOn = mean(cell2mat(firingRatesPrefDir));

    %Finding + caclulating average max FR for ND off 
    nullDirIndices = find(stimDir == absNullDirOn);
    for y = 1:numel(nullDirIndices)
        index = nullDirIndices(y);
        maxFR = maxFRon(index);
        firingRatesNullDir{y} = maxFR;
    end
%     disp(firingRatesNullDir);
    avgFRnullOn = mean(cell2mat(firingRatesNullDir));

    %on avg FR DSI

    DSIonFR = (avgFRprefOn - avgFRnullOn)/(avgFRprefOn + avgFRnullOn);

% Initialize storage for preferred and null direction firing rates
onFiringRatesPD = [];
onFiringRatesND = [];

% Finding and calculating average FR for preferred direction (PD) ON
prefDirIndices = find(stimDir == absPrefDirOn);
for y = 1:numel(prefDirIndices)
    index = prefDirIndices(y);
    onFiringRatesPD = [onFiringRatesPD; firingRatesOnWindow(index, :)];
end


% Finding and calculating average FR for null direction (ND) ON
nullDirIndices = find(stimDir == absNullDirOn);
for y = 1:numel(nullDirIndices)
    index = nullDirIndices(y);
    onFiringRatesND = [onFiringRatesND; firingRatesOnWindow(index, :)];
end

avgPDfiringRateON = mean(onFiringRatesPD, 1);
avgNDfiringRateON = mean(onFiringRatesND, 1);

%For OFF

 for x = 1:24

   spCts = spCtsOFF(x,:); 
   spCts = ones(1,spCts)'; 
   spTms = spTmsOFF{x,:}';
   spikenTimes = [spCts spTms];

   %define parameters for sliding window 
   windowSize = 0.1; %100ms sliding window
   totalTime = 1.5; %total duration of spike phase 
   numWindows = floor(totalTime / windowSize); 

   %Variables to store firing rates and window start times
   firingRates = zeros(numWindows,1); 
   windowStartTimes = (0:numWindows-1) * windowSize;
   allFiringRates = zeros(15, numWindows);

   %Sliding the window 
   for j = 1:15

    startTime = (j - 1) * windowSize; 
    endTime = startTime + windowSize; 

    %count number of spikes in current window 
    spikesInWindow = sum(spikenTimes(:,2) >= startTime & spikenTimes(:,2) < endTime);
    %firing rate for current window
    firingRates(j) = spikesInWindow /windowSize; 
   end 
   
    %store all firing rates for each window 
    allFiringRates = max(allFiringRates, firingRates); 
    allFiringRates = allFiringRates(:,1);
    firingRatesOffWindow(x,:) = firingRates; 

    maxFiringRate = max(allFiringRates);
    maxFRoff(x) = maxFiringRate; 
    
 end 

%  disp(maxFRoff); 

    %Finding + caclulating average max FR for PD off 
    prefDirIndices = find(stimDir == absPrefDirOff);
    for y = 1:numel(prefDirIndices)
        index = prefDirIndices(y);
        maxFR = maxFRoff(index);
        firingRatesPrefDir{y} = maxFR;
    end

    avgFRprefOff = mean(cell2mat(firingRatesPrefDir));

    %Finding + caclulating average max FR for ND off 
    nullDirIndices = find(stimDir == absNullDirOff);
    for y = 1:numel(nullDirIndices)
        index = nullDirIndices(y);
        maxFR = maxFRoff(index);
        firingRatesNullDir{y} = maxFR;
    end

    avgFRnullOff = mean(cell2mat(firingRatesNullDir));

    DSIoffFR = (avgFRprefOff - avgFRnullOff)/(avgFRprefOff + avgFRnullOff);

% Initialize storage for preferred and null direction firing rates
offFiringRatesPD = [];
offFiringRatesND = [];

% Finding and calculating average FR for preferred direction (PD) Off
prefDirIndices = find(stimDir == absPrefDirOff);
for y = 1:numel(prefDirIndices)
    index = prefDirIndices(y);
    offFiringRatesPD = [offFiringRatesPD; firingRatesOffWindow(index, :)];
end


% Finding and calculating average FR for null direction (ND) OFF
nullDirIndices = find(stimDir == absNullDirOff);
for y = 1:numel(nullDirIndices)
    index = nullDirIndices(y);
    offFiringRatesND = [offFiringRatesND; firingRatesOffWindow(index, :)];
end

avgPDfiringRateOFF = mean(offFiringRatesPD, 1);
avgNDfiringRateOFF= mean(offFiringRatesND, 1);


%% Functions
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
    dTable.barSpeed(i) = expTable.barSpeed(i);
    dTable.StimDir(i,:) = stimDir';
    dTable.SpikeprefDirTrace{i} = SpikeprefDirTrace;
    dTable.SpikenullDirTrace{i} = SpikenullDirTrace;
    dTable.DSIon(i) = DSIon;
    dTable.DSIoff(i) = DSIoff;
    dTable.normVSon(i) = normVSon; 
    dTable.normVSoff(i) = normVSoff; 
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
    dTable.maxFRon{i} = maxFRon; 
    dTable.maxFRoff{i} = maxFRoff; 
    dTable.allMaxFR{i} = allMaxFR; 
    dTable.avgFRprefOn(i) = avgFRprefOn; 
    dTable.avgFRprefOff(i) = avgFRprefOff; 
    dTable.avgFRnullOn(i) = avgFRnullOn; 
    dTable.avgFRnullOff(i) = avgFRnullOff; 
    dTable.firingRatesAllWindows{i} = firingRatesAllWindows;
    dTable.avgPDfiringRate{i} = avgPDfiringRate;
    dTable.avgNDfiringRate{i} = avgNDfiringRate;
    dTable.avgPDfiringRateON{i} = avgPDfiringRateON; 
    dTable.avgNDfiringRateON{i} = avgNDfiringRateON; 
    dTable.avgPDfiringRateOFF{i} = avgPDfiringRateOFF; 
    dTable.avgNDfiringRateOFF{i} = avgNDfiringRateOFF;
    dTable.DSIonFR(i) =DSIonFR;
    dTable.DSIoffFR(i) =DSIoffFR;
    dTable.PDdiff(i) = PDdiff; 

  
end


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
    % axis tight
    ylim([-100 100])
    
end
title(titleStr);
end

% Function to process each phase
function [spTms, spCts] = process_phase(data, thresh, si, manual)
    % Spike counts for given phase
    [nPts, nTrls] = size(data);
    dVector = reshape(data, numel(data), 1);
    dFilter = filterData(dVector, si);
    thresh = abs(thresh);
    spTmsRaw = getSpikeTimesCurrentClamp(dFilter, si, 'thresh', thresh, 'manual', manual);
    secTrial = nPts * si * 1e-6; % Convert to seconds
    spTms = cell(nTrls, 1);
    spCts = NaN(nTrls, 1);

    for f = 1:nTrls
        spIndx = spTmsRaw > (f-1)*secTrial & spTmsRaw < f*secTrial;
        spikes = spTmsRaw(spIndx) - (f-1)*secTrial;
        spTms{f} = spikes;
        spCts(f) = numel(spikes);
    end
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

dsgcDir = '/Users/karinabistrong/Dropbox/Mac/Desktop/DATA/Ephys';
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
preSet_rlim = [0 100];
prefSpikes = ctMean(prefIndx);
nullSpikes = ctMean(nullIndx);
titleStr = sprintf('Pref Dir %3.1f DSI %4.2f Vec Length %4.2f',prefDir,DSI,vecLength);
acqInfo = sprintf('Neuron %g Date %g Hold %g Drug %s',cellNumb,expDate,holdVol,drugCond);
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
    
    title(titleStr,acqInfo);
    rlim(preSet_rlim);
%     if rFlag
%         pAx.RLim = rBounds;
%     end
    
% else
    hF = [];
end
end


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
        'YLim',[-100 100],'XLim',[tSec(1) tSec(end)],'NextPlot','replacechildren');
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
    ylim([-100 100])
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