%% Build the neuron table
%
% This code needs an excel spreadsheet with metadata in it.
% Right now, it will read 'EphysExperimentSheet.xlsx'
% Every row in that excel spreadsheet is a single recording
% The buildDTable code will run all the recording
%
% The end result of this code is that it saves 'dTable.mat'
% This is the table where each row is a recording
% Other codes will load this table and run analysis on it

%% Load the table
data_guide_name = 'EphysExperimentSheet.xlsx';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'ExperimentDate','path','FileName','Sex','Age', 'Genotype', 'GFPLabel',...
    'Eye', 'Quadrant','RecordingType','Stim','StimFile','StimTimes','CellNumber','Drugs','barLength','barSpeed','numberOfRecs'},'string');

expTable = readtable(data_guide_name, opts');

[num_files, dummy] = size(expTable);

%% Intialize the neuron table (called neuronTable)

totalRecs = height(expTable);
neuronCounter = 1; %This value will keep track of which neuron I am on

dTable = table('Size', [totalRecs 16], 'VariableTypes', {'double','double','string','string','string',...
    'string','string','string','string','double','string','string','string','double','string','double'});

dTable.Properties.VariableNames = {'expDate', 'aquisitionNumber', 'sex','age', 'genotype',...
    'GFPLabel', 'eye', 'quad','recType','holdVol','stim','stimFile','stimTimes','cellNumb','drugs','holdCur'};

dTable.StimDir = NaN(totalRecs,24);
dTable.pCurOn = NaN(totalRecs,24); 
dTable.pCurOff = NaN(totalRecs,24); 

%% For loop for going through each recording to populate table

for i = 1:2  % this will run through all files
% for i = 1:2 % For testing code, num_files for all files
    
    cd (expTable.path(i))

% making columns for these variables 
    stimDir = nan(1,24);
    pCurON = nan(1,24);
    pCurOFF = nan(1,24);
    DSIon = nan(1);
    DSIoff = nan(1);
    prefDirOn = nan(1);
    prefDirOff = nan(1);
    
     %for bars
    if expTable.Stim(i) == "bars"
       %load abf file
       [d,si]=abfload(expTable.FileName(i));
       %load raw currents
       if expTable.AquisitionNumber(i) < 10 
        quickLoad(append("00",num2str(expTable.AquisitionNumber(i))));
       end 
       if expTable.AquisitionNumber(i) > 9 
        quickLoad(append("0",num2str(expTable.AquisitionNumber(i))));
       end 
       %load stim info to get stim times
       if expTable.AquisitionNumber(i) < 10 
       quickStim(append("00",num2str(expTable.AquisitionNumber(i))));
       end 
       if expTable.AquisitionNumber(i) > 9 
         quickStim(append("0",num2str(expTable.AquisitionNumber(i))));  
       end 
       %select just the stimulus direction
       stimDir = stimInfo(1:24,1);

        %normalize baseline +lowpass filter 
        if expTable.holdCur(i) ~= NaN
            d = d - expTable.holdCur(i); 
            dt = 1e-4;
            d = lowpass(d,30,1/dt,'ImpulseResponse','iir','Steepness',.8); % lowpass filter
        end   

        %get peak amplitudes for inhibitory current 
        if expTable.holdVol(i) > 0
            CurON = d(10000:15000,1:24);
            pCurON = max(CurON);
            CurOFF = d(20000:25000,1:24);
            pCurOFF = max(CurOFF);
        end
        
        %peak amplitude for excitatory currents
        if expTable.holdVol(i) < 0
            CurON = d(10000:15000,1:24);
            pCurON = min(CurON);
            CurOFF = d(20000:25000,1:24);
            pCurOFF = min(CurOFF);
        end 
        
        %sort based on amplitude and direction for both ON and OFF 
        pSortON = quickSort(pCurON,stimDir);
        pSortOFF = quickSort(pCurOFF,stimDir);
        
        %plot tuning curve and get DSI and pref dir for ON and OFF 
        [prefDir, DSI] = dirTuning(pSortON,stimDir);
        DSIon = DSI;
        prefDirOn =prefDir;
        [prefDir, DSI,vecLength] = dirTuning(pSortOFF,stimDir);
        DSIoff = DSI;
        prefDirOff = prefDir;
        
    end 
    
    % Load neuron info in the neuronTable
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
    dTable.stimTimes(i) = expTable.StimTimes(i);
    dTable.stimFile(i) = expTable.StimFile(i);
    dTable.cellNumb(i) = expTable.CellNumber(i);
    dTable.drugs(i) = expTable.Drugs(i);
    dTable.holdCur(i) = expTable.holdCur(i);
    dTable.quad(i) = expTable.Quadrant(i);
    dTable.StimDir(i,:) = stimDir';
    dTable.pCurOn(i,:) = pCurON';
    dTable.pCurOff(i,:) = pCurOFF';
    dTable.DSIon(i) = DSIon;
    dTable.DSIoff(i) = DSIoff;
    dTable.prefDirOn(i) = prefDirOn;
    dTable.prefDirOff(i) = prefDirOff;
end 

%% Save table 

save('/Users/kbistrong/Dropbox/My Mac (karina’s MacBook Pro)/Desktop/DATA/Ephys/compiledResults/dTable','dTable');


%% Functions 

function hF = quickLoad(abfStim,dirName,hideOutput)
%Simple function to load an abf file and compress singleton dimensions.
%If provided only an abf file stim number, loads from current directory.
%Will automatically display raw data unless suppressed. Requires abfload.
%Inputs
%   abfStim         string, last 3 digits of abf file to be loaded
%   dirName         string, directory to load from. Assumes current
%   directory if not provided. Function interprets beginning of abf file
%   name from this string, so directory name must be date of recording in
%   Clampex's date format.
%   hideOutput      true/false, flag to hide raw data plots.
%Outputs
%   hF              figure handle of raw data trace
%   d               matrix, number of pts x number of trials (assigned in)
%   si              int, sampling interval in microseconds (assigned in)
% MTS 5/15/2016

if nargin < 3
    %Unless requested, don't display raw data traces
    hideOutput = false;
end

if nargin < 2 || isempty(dirName) 
    %if no input argument, assume current directory.
    [~,dirName] = fileparts(pwd);
    fprintf('Loading from current folder; %s \n',dirName);
else
    dirName = char(dirName); %ensure indexable characters, not a string
end

assert(size(abfStim,1) == 1,'abfStim must be a single string, not an array of strings.');

%account for clampex enforcing 5 digit dates
if strcmp(dirName(3:4),'10')
    abfDate = [dirName(1:2) 'o' dirName(5:6)];
elseif strcmp(dirName(3:4),'11')
    abfDate = [dirName(1:2) 'n' dirName(5:6)];
elseif strcmp(dirName(3:4),'12')
    abfDate = [dirName(1:2) 'd' dirName(5:6)];
elseif strcmp(dirName(3),'0') %account for clampex's peculiar naming conventions
    abfDate = [dirName(1:2) dirName(4:6)];
else
    abfDate = dirName(1:5);
end

abfName = sprintf('%s%s.abf',abfDate,abfStim);

%following line assumes Karina's computer and system architecture
dsgcDir = '/Users/kbistrong/Dropbox/My Mac (karina’s MacBook Pro)/Desktop/DATA/Ephys';
searchDirName = sprintf('%s*',dirName); %find directories that match input date

oldDir = cd(dsgcDir);
newDir = dir(searchDirName);
cd(newDir.name);

[d,si] = abfload(abfName);
d = squeeze(d);


%%%%uncomment below to plot raw traces 
if ~hideOutput
    size(d)
    hF = plotRawData(d,si);
else
    hF = [];
end

cd(oldDir);

%dispense variables to caller function, so don't have to declare outputs
%each time I'm calling quickLoad
assignin('caller','d',d);
assignin('caller','si',si);

end

function hF = plotRawData(d,si)
%plots unprocessed data if called, outputs figure handle.

[nPts,trials] = size(d);
dt = si*1e-6;
t = 0:dt:(nPts - 1)*dt;
L = floor(sqrt(trials));
W = ceil(trials / L);
minD = min(d(:));
maxD = max(d(:));
hF = figure;
for i = 1:trials
    subplot(L,W,i)
    plot(t,d(:,i),'r') 
    axis tight
    ylim([minD maxD])
end

end

function quickStim(stimNum,stimDate)
% Function to load stimXXX.txt files for a given experiment.

if nargin < 2 || isempty(stimDate) 
    %if no input argument, assume current directory.
    [~,stimDate] = fileparts(pwd);
    fprintf('Loading from current folder; %s \n',stimDate);
else
    stimDate = char(stimDate); %ensure indexable characters, not a string
end

dsgcDir = '/Users/kbistrong/Dropbox/My Mac (karina’s MacBook Pro)/Desktop/DATA/Ephys';
searchDirName = sprintf('%s*',stimDate); %find directories that match input date

oldDir = cd(dsgcDir);
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
nuX = numel(uX);
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

function [prefDir, DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(ctSort,stimDirs)
% add HF back in and uncomment below in order to have polar plots
% generated
if nargin<5 || isempty(rBounds)
    rFlag = false;
else
    rFlag = true;
end

if nargin<4 || isempty(makeFig)
    makeFig = 1;
end
if nargin<3 || isempty(showLess)
    showLess = 0;
end

ctMean = mean(ctSort,2);
% ctSortPlot = [ctSort; ctSort(1,:)];
% ctMeanPlot = [ctMean; ctMean(1,:)];
% 
% nReps = size(ctSort,2);

uDirs = unique(stimDirs); %potentially problematic way of doing things;
% if stimDirs is shifted in a non-uniform way (e.g. resets to 0 after passing 360)
% then ctSort will now be incorrectly indexed
uDirs = deg2rad(uDirs);

[x,y] = pol2cart(uDirs, ctMean / sum(ctMean));
prefDir = atan2d(sum(y),sum(x));
if prefDir < 0
    prefDir = prefDir + 360;
end
vecLength = sqrt(sum(x)^2 + sum(y)^2);

nDirs = length(uDirs);
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

prefSpikes = ctMean(prefIndx);
nullSpikes = ctMean(nullIndx);
% % uncomment everything below to get polar plots of each recording 
% titleStr = sprintf('Pref Dir %3.1f DSI %4.2f Vec Length %4.2f',prefDir,DSI,vecLength);
% 
% uDirs = [uDirs; uDirs(1)];
% uDirsPlot = repmat(uDirs,1,nReps);
% 
% 
% if prefDir < 0
%     prefDir = prefDir + 360;
% end
% 
% %% Plot figure
% 
%  if ~showLess
%     if makeFig
%         hF = figure;
%     else
%         hF = [];
%     end
%     
%     %Set radial bounds
%     pAx = polaraxes();
%     
%     plotChild = polarplot(uDirsPlot,ctSortPlot);
%     set(plotChild(1:nReps),'LineWidth',1);
%     
%     hold on
%     
%     plotChild = polarplot(uDirs,ctMeanPlot,'k');
%     set(plotChild(1),'LineWidth',2)
%     
%     %plotChild = compass(sum(x)*max(ctMean),sum(y)*max(ctMean),'k'); %think about this
%     plotChild = polarplot(deg2rad([prefDir, prefDir]), [0 max(ctMean)*vecLength],'k');
%     set(plotChild(1),'LineWidth',1.5)
%     
%     title(titleStr);
%     if rFlag
%         pAx.RLim = rBounds;
%     end
%     
% else
%     hF = [];
% 
%  end 
end 
