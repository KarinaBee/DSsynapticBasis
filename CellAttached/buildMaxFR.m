%% Build the neuron table FOR CELL ATTACHED DATA 

%This code just loads the data into cTable WITH PLOTS. See commented sections below 
% if on the macbook- REMEMBER to change the directory to kbistrong instead of karina bistrong & in the master spreadsheet 

% This code needs an excel spreadsheet with metadata in it.

% Every row in that excel spreadsheet is a single recording
% The buildTable code will run all the recording
%
% The end result of this code is that it saves 'cTable.mat'
% This is the table where each row is a recording
% Other codes will load this table and run analysis on it

%% Load the table
data_guide_name = 'B2vertCellAttached.xlsx';

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
    % cd (expTable.path(i))
    cd (expTable.macPath(i))

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
       Eye = expTable.Eye(i);

       [newInd] = eyeCalSOS(stimDir,Eye); 

       stimDirOld = stimDir;
       stimDir(:,1) = newInd; 

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

        % Parameters
        thresh = -7; % Default threshold value for spike detection
        manual = 0; % Manual flag

        % Process ON phase
        [spTmsON, spCtsON] = process_phase(dON, thresh, si, manual);
        assignin('caller', 'spTmsON', spTmsON);
        assignin('caller', 'spCtsON', spCtsON);

        % Process OFF phase
        [spTmsOFF, spCtsOFF] = process_phase(dOFF, thresh, si, manual);
        assignin('caller', 'spTmsOFF', spTmsOFF);
        assignin('caller', 'spCtsOFF', spCtsOFF);

        % Define parameters for sliding window 
        windowSize = 0.1; % 100ms window
        totalTime = 1.5;  % total duration of ON phase
        numWindows = floor(totalTime / windowSize); 
        
        % Preallocate for efficiency (keeps structure intact)
        firingRatesOnWindow = zeros(24, numWindows); 
        
        %% FOR ON 
        
        % Loop through all 24 trials
        for x = 1:24
            spTms = spTmsON{x}; % Get spike times for this trial

        % Compute firing rates using the sliding window
        for j = 1:numWindows
            startTime = (j - 1) * windowSize; 
            endTime = startTime + windowSize; 

        % Count spikes in window (empty trials naturally remain 0)
        firingRatesOnWindow(x, j) = sum(spTms >= startTime & spTms < endTime) / windowSize; % Hz
        end 

        end

        maxFRon = max(firingRatesOnWindow, [], 2); % Max firing rate per trial (24x1)
        pSortON = quickSort(maxFRon,stimDir);
        dirList = [0; 45; 90; 135; 180; 225; 270; 315];
        [hF, prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSortON,stimDir,cellNumb,expDate,holdVol,drugCond);

        prefDirOn = prefDir; 
        DSIon =DSI; 
        prefFRon = prefSpikes; 
        nullFRon = nullSpikes; 

        if prefDirOn > 337.5
            prefDirOn = 0;
        end 

        [temp, prefInd] = min(abs(dirList - prefDir));
        absPrefDirOn = dirList(prefInd); 
        absNullDirOn = absPrefDirOn - 180;
        if absNullDirOn < 0 
            absNullDirOn = absNullDirOn + 360;
        end

%% FOR OFF 
        
% Loop through all 24 trials
for x = 1:24
    spTms = spTmsOFF{x}; % Get spike times for this trial

    % Compute firing rates using the sliding window
    for j = 1:numWindows
        startTime = (j - 1) * windowSize; 
        endTime = startTime + windowSize; 

        % Count spikes in window (empty trials naturally remain 0)
        firingRatesOffWindow(x, j) = sum(spTms >= startTime & spTms < endTime) / windowSize; % Hz
    end 
end

maxFRoff = max(firingRatesOffWindow, [], 2); % Max firing rate per trial (24x1)
pSortOFF = quickSort(maxFRoff,stimDir);
dirList = [0; 45; 90; 135; 180; 225; 270; 315];
[hF, prefDir, DSI, vecLength, prefSpikes, nullSpikes] = dirTuning(pSortOFF, stimDir, cellNumb, expDate, holdVol, drugCond);

prefDirOff = prefDir; 
DSIoff = DSI; 
prefFRoff = prefSpikes; 
nullFRoff = nullSpikes; 

if prefDirOff > 337.5
    prefDirOff = 0;
end 

[temp, prefInd] = min(abs(dirList - prefDir));
absPrefDirOff = dirList(prefInd); 
absNullDirOff = absPrefDirOff - 180;
if absNullDirOff < 0 
    absNullDirOff = absNullDirOff + 360;
end

%% Saving to dTable 
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
    dTable.DSIon(i) = DSIon;
    dTable.DSIoff(i) = DSIoff;
    dTable.absPrefDirOn(i) = absPrefDirOn;
    dTable.absPrefDirOff(i) = absPrefDirOff;
    dTable.absNullDirOn(i) = absNullDirOn;
    dTable.absNullDirOff(i) = absNullDirOff;
    dTable.prefFRon(i) = prefFRon; 
    dTable.nullFRon(i) = nullFRon; 
    dTable.prefFRoff(i) = prefFRoff; 
    dTable.nullFRoff(i) = nullFRoff;
end 

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
    % axis tight
    ylim([-100 100])
    
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
