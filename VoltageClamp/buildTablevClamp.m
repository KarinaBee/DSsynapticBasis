%% Build the neuron table

%This code loads the data into dTable WITH PLOTS. See commented sections below 

% This code needs an excel spreadsheet with metadata in it.
% Right now, it will read 'EphysExperimentSheet.xlsx'
% Every row in that excel spreadsheet is a single recording
% The buildTable code will run all the recording
%
% The end result of this code is that it saves 'dTable.mat'
% This is the table where each row is a recording
% Other codes will load this table and run analysis on it
%functions used are at the bottom of teh script

%% Load the table
data_guide_name = '';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);


% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'ExperimentDate','macPath','path','FileName','Sex','Age', 'Genotype', 'GFPLabel',...
    'Eye', 'Quadrant','RecordingType','Stim','StimFile','StimTimes','CellNumber','Drugs','barLength','barSpeed','numberOfRecs'},'string');

expTable = readtable(data_guide_name, opts');

[num_files, dummy] = size(expTable); %determines number of rows (iterations for loop will go through)
 
%creating empty lists for inhibitory PD for excitatory PD to use
inhDirOn = nan(num_files, 1); 
resultOn = length(inhDirOn); 
inhDirOff = nan(num_files, 1);
resultOff = length(inhDirOff); 

%% Intialize the neuron table (called dTable)

totalRecs = height(expTable);
neuronCounter = 1; %This value will keep track of which neuron I am on

dTable = table('Size', [totalRecs 16], 'VariableTypes', {'double','double','string','string','string',...
    'string','string','string','string','string','string','string','string','double','string','double'});

dTable.Properties.VariableNames = {'expDate', 'aquisitionNumber', 'sex','age', 'genotype',...
    'GFPLabel', 'eye', 'quad','recType','holdVol','stim','stimFile','stimTimes','cellNumb','drugs','holdCur'};

dTable.StimDir = NaN(totalRecs,24);
dTable.pCurOn = NaN(totalRecs,24); 
dTable.pCurOff = NaN(totalRecs,24); 
dTable.pTimeOn = NaN(totalRecs,24); 
dTable.pTimeOff = NaN(totalRecs,24); 
dTable.avgPsortON = NaN(totalRecs,8); 
dTable.avgPsortOFF = NaN(totalRecs,8);
dTable.prefDirTrace = cell(totalRecs, 1); % Initialize as a cell array
dTable.nullDirTrace = cell(totalRecs, 1); % Initialize as a cell array

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

% making columns for these variables 
    stimDir = nan(1,24);
    pCurON = nan(1,24);
    pCurOFF = nan(1,24);
    pTimeOn = nan(1,24);
    pTimeOff = nan(1,24);
    DSIon = nan(1);
    DSIoff = nan(1);
    normVSon = nan(1);
    normVSoff = nan(1); 
    avgPsortOFF = nan(8,1);
    prefDirOn = nan(1);
    prefDirOff = nan(1);
    absPrefDirOn = nan(1);
    absPrefDirOff = nan(1);
    absNullDirOn = nan(1);
    absNullDirOff = nan(1);
    pCurPrefON = nan(1);
    pCurPrefOFF = nan(1);
    pCurNullON = nan(1);
    pCurNullOFF = nan(1);
    contrast = nan(1);
    Rin = nan(1); 
    eAsymON = nan(1);
    eAsymOFF = nan(1);
    iAsymON = nan(1);
    iAsymOFF = nan(1);
    
    %for bars
    if expTable.Stim(i) == "bars" 

     %load abf file
       [d,si]=abfload(expTable.FileName(i));
  
        %normalize baseline + lowpass filter for negative holdCur 
            x = mean(d(1:3000,:,:)); %averaging across the first .3 seconds 
            d_corr = d(:,:,:) - x;
            dt = 1e-4; %sampling rate 10Hz 
            d =squeeze(d_corr);
            d  = lowpass(d,30,1/dt,'ImpulseResponse','iir','Steepness',.8); %lowpass filter 
           
%plots filtered raw data, comment to suppress raw data plots 
       hF = rawTest(d,si,cellNumb,GFPlabel,expDate,holdVol);
    

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

%% Get peak amplitudes and times

% for inhibitory current 
        if expTable.holdVol(i) > 0
            CurON = d(10000:22000,1:24);
            pCurON = max(CurON); 
            [pCurON, pTimeON] = max(CurON); %times of peak current
            pTimeON = pTimeON + 10000;
            CurOFF = d(pTimeON+5000:40000,1:24);
            pCurOFF = max(CurOFF);
            [pCurOFF, pTimeOFF] = max(CurOFF); %times of peak current
            pTimeOFF = pTimeOFF + (pTimeON+5000);
        end
        
%for excitatory currents
        if expTable.holdVol(i) < 0
            CurON = d(15000:25000,1:24);
            pCurON = min(CurON);
            [pCurON, pTimeON] = min(CurON); %times of peak current
            pTimeON = pTimeON + 15000;
            CurOFF = d(pTimeON+5000:40000,1:24);
            pCurOFF = min(CurOFF);
            [pCurOFF, pTimeOFF] = min(CurOFF); %times of peak current
            pTimeOFF = pTimeOFF + (pTimeON+5000); 
        end 

%% Sort data for ON and OFF based on direction & plot them
        pSortON = quickSort(pCurON,stimDir);
        pSortOFF = quickSort(pCurOFF,stimDir);

        nDirectionsON = height(pSortON);
        avgPsortON = [];
       
        for bb = 1:nDirectionsON
            rowAverage = mean(pSortON(bb,:));
            avgPsortON = [avgPsortON; rowAverage];
        end 

        nDirectionsOFF = height(pSortOFF);
        avgPsortOFF = [];
        for bb = 1:nDirectionsOFF
            rowAverage = mean(pSortOFF(bb,:));
            avgPsortOFF = [avgPsortOFF; rowAverage];
        end 

        
        [avgTrace,prefDirTrace,nullDirTrace] = plotDirTraces(d,stimDir,1,10000,absNullDirOn,absPrefDirOn);

        prefDirTrace = mean(prefDirTrace,2);
        nullDirTrace = mean(nullDirTrace,2);
        pdt = size(prefDirTrace,1);
        ndt = size(nullDirTrace,1);
       
        dTable.prefDirTrace = NaN(totalRecs,pdt);
        dTable.nullDirTrace = NaN(totalRecs,ndt);

        %normalized VS for inhibition
        if expTable.holdVol(i) > 0 
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
       end 

        %normalized VS for inhibition
        if expTable.holdVol(i) < 0 
            %for ON
            pMeanON = mean(pSortON,2); %average response for each direction 
            pMaxON = min(pMeanON); %finds max response
            pNormON = pMeanON/pMaxON; %normalizes all resp to max resp
            %for OFF
            pMeanOFF = mean(pSortOFF,2); 
            pMaxOFF = min(pMeanOFF); 
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
       end
%% get DSI, vec length and pref dir for ON and OFF - uncomment in function to plot
        [hF,prefDir,DSI,vecLength] = dirTuning(pSortON,stimDir,cellNumb,expDate,holdVol,drugCond);
        DSIon = DSI;
        prefDirOn =prefDir;
        vecSumOn = vecLength;

        %calculating and storing inh PD 
        for jj = 1:resultOn 
            inhDirOn(i) = prefDirOn; 
        end 
        %if excitatory, replacing PD with indexed inh PD 
        if expTable.holdVol(i) < 0 
            prefDirOn = inhDirOn(i-1);
        end 
       
        [hF,prefDir,DSI,vecLength] = dirTuning(pSortOFF,stimDir,cellNumb,expDate,holdVol,drugCond);
        DSIoff = DSI;
        vecSumOff = vecLength; 
        prefDirOff = prefDir;

        %same as above but for Off component 
        for kk = 1:resultOff
            inhDirOff(i) = prefDirOff; 
        end 

        if expTable.holdVol(i) < 0 
            prefDirOff = inhDirOff(i-1);
        end 

%% get pref and null dir from shown direction list 
    % FOR ON
        if prefDirOn > 337.5
            prefDirOn = 0;
        end 
        [temp, prefInd] = min(abs(stimDir - prefDirOn));
        absPrefDirOn = stimDir(prefInd); 
        absNullDirOn = absPrefDirOn - 180;
        if absNullDirOn < 0 
            absNullDirOn = absNullDirOn + 360;
        end 

    % FOR OFF
        if prefDirOff > 337.5
            prefDirOff = 0;
        end 
        [temp, prefInd] = min(abs(stimDir - prefDirOff));
        absPrefDirOff = stimDir(prefInd); 
        absNullDirOff = absPrefDirOff - 180; 
        if absNullDirOff < 0 
            absNullDirOff = absNullDirOff + 360;
        end 

        [avgTrace,prefDirTrace,nullDirTrace] = plotDirTraces(d,stimDir,1,10000,absNullDirOn,absPrefDirOn);
       % 
       % [avgTrace,prefDirTrace,nullDirTrace] = plotDirTraces(dFiltered,stimDir,1,10000,absNullDirOn,absPrefDirOn);

        prefDirTrace = mean(prefDirTrace,2);
        nullDirTrace = mean(nullDirTrace,2);
%% get average peak response for pref and null dir
    x = stimDir; 
    y = pCurON'; 
    z = pCurOFF';

    %for pref ON
    j = find(x==absPrefDirOn);
    k = y(j); 
    pCurPrefON = mean(k); 
    %for null ON
    l = find(x==absNullDirOn);
    m =y(l);
    pCurNullON = mean(m);
    %for pref OFF
    n = find(x==absPrefDirOff);
    o = z(n);
    pCurPrefOFF = mean(o);
    %for null OFF
    p = find(x==absNullDirOff);
    q = z(p);
    pCurNullOFF = mean(q); 
 

%get average peak time for preferred and null directions 
    r = pTimeON'; 
    s = pTimeOFF';

 %for pref ON
    t = find(x==absPrefDirOn);
    u = r(t); 
    pTimePrefON = mean(u); 
    %for null ON
    v = find(x==absNullDirOn);
    w = r(v);
    pTimeNullON = mean(w);
    %for pref OFF
    a = find(x==absPrefDirOff);
    b = s(a);
    pTimePrefOFF = mean(b);
    %for null OFF
    c = find(x==absNullDirOff);
    e = s(c);
    pTimeNullOFF = mean(e);


    %get asymmetric inhibitory component 
    if expTable.holdVol(i) < 0 
       eAsymON = abs(pCurNullON - pCurPrefON); 
       eAsymOFF = abs(pCurNullOFF - pCurPrefOFF); 
    end

    if expTable.holdVol(i) > 0 
        iAsymON = abs(pCurPrefON - pCurNullON); 
        iAsymOFF = abs(pCurPrefOFF - pCurNullOFF);
    end
    
%%  Load neuron info in the neuronTable
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
    dTable.pTimeOn(i,:)= pTimeON';
    dTable.pTimeOff(i,:)= pTimeOFF';
    dTable.avgPsortON(i,:) = avgPsortON';
    dTable.avgPsortOFF(i,:) = avgPsortOFF';
    dTable.DSIon(i) = DSIon;
    dTable.DSIoff(i) = DSIoff;
    dTable.normVSon(i) = normVSon; 
    dTable.normVSoff(i) = normVSoff; 
    dTable.vecSumOn(i) = vecSumOn;
    dTable.vecSumOff(i) =vecSumOff; 
    dTable.prefDirOn(i) = prefDirOn;
    dTable.prefDirOff(i) = prefDirOff;
    dTable.absPrefDirOn(i) = absPrefDirOn;
    dTable.absPrefDirOff(i) = absPrefDirOff;
    dTable.absNullDirOn(i) = absNullDirOn;
    dTable.absNullDirOff(i) = absNullDirOff;
    dTable.pCurPrefON(i) = pCurPrefON;
    dTable.pCurPrefOFF(i) = pCurPrefOFF;
    dTable.pCurNullON(i) = pCurNullON;
    dTable.pCurNullOFF(i) = pCurNullOFF;
     dTable.eAsymON(i) = eAsymON; 
    dTable.eAsymOFF(i) = eAsymOFF;
    dTable.iAsymON(i) = iAsymON; 
    dTable.iAsymOFF(i) = iAsymOFF;
    dTable.pTimePrefOFF(i) = pTimePrefOFF;
    dTable.pTimePrefON(i) = pTimePrefON;
    dTable.pTimeNullON(i) = pTimeNullON;
    dTable.pTimeNullOFF(i) = pTimeNullOFF;
    dTable.contrast(i) = expTable.contrast(i);
 
    end
    end


% cd 'C:\Users\Karina Bistrong\Dropbox\Mac\Desktop\DATA\Ephys\compiledResults\WholeCell'
% save('dTable.mat','dTable');
%% Functions 

function hF = rawTest(d,si,cellNumb,GFPlabel,expDate,holdVol)
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

titleStr = sprintf('Neuron %g GFP %g Date %g Hold %g',cellNumb,GFPlabel,expDate,holdVol);

for i = 1:trials
    subplot(L,W,i)
    plot(t,d(:,i),'r') 
    axis tight
    ylim([minD 2500])
    
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

function [avgTrace,prefDirTrace, nullDirTrace] = plotDirTraces(stimDF,stimDir,ROI,Fs,absNullDirOn,absPrefDirOn)

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
        'YLim',[mindF 2000],'XLim',[tSec(1) tSec(end)],'NextPlot','replacechildren');
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
    ylim([mindF 2000])
    xlim([tSec(1) tSec(end)])
    title([ num2str(uDirs(j))])
    txtPos = [(xStart + (j-1)*xEach) yEnd xEach yLabelWidth];
    txtStr = num2str(uDirs(j));
    txt = uicontrol(hF,'Style','text','Units','Normalized','Position',txtPos,...
        'String',txtStr,'FontSize',12);
end   


% Extract the subplot for the preferred direction "0"
prefDirIndex = find(uDirs == absPrefDirOn); % Find the index of the preferred direction "0"
prefDirTrace = dFSort(:, 1:nReps, prefDirIndex); % Extract the data for the preferred direction
% Extract the subplot for the null direction "180"
nullDirIndex = find(uDirs == absNullDirOn); % Find the index of the null direction "180"
nullDirTrace = dFSort(:, 1:nReps, nullDirIndex); % Extract the data for the null direction

end

function [hF,prefDir,DSI,vecLength,prefSpikes,nullSpikes] = dirTuning(pSort,stimDir,cellNumb,expDate,holdVol,drugCond)

% if nargin<5 || isempty(rBounds)
%     rFlag = false;
% else
%     rFlag = true;
% end
% 
% if nargin<4 || isempty(makeFig)
%     makeFig = 1;
% end
% if nargin<3 || isempty(showLess)
%     showLess = 0;
% end

ctMean = mean(pSort,2); %takes the mean of each row in the sorted ON or OFF responses = average response for each direction 
ctSortPlot = [pSort; pSort(1,:)]; %produces a matrix of the sorted responses ie is just pSortON or pSort OFF
ctMeanPlot = [ctMean; ctMean(1,:)]; %matrix of the means 

nReps = size(pSort,2); %defining number of repitions of directions based on number of columns in pSort 

uDirs = unique(stimDir);% list of directions %potentially problematic way of doing things;
% if stimDirs is shifted in a non-uniform way (e.g. resets to 0 after passing 360)
% then ctSort will now be incorrectly indexed
uDirs = deg2rad(uDirs); %converts list of directions to radians 
% [x,y] = pol2cart(uDirs, ctMean);
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
preSet_rlim = [0 1500];
prefSpikes = ctMean(prefIndx);
nullSpikes = ctMean(nullIndx);
titleStr = sprintf('Pref Dir %3.1f DSI %4.2f Vec Length %4.2f',prefDir,DSI,vecLength);
acqInfo = sprintf('Neuron %g Date %g Hold %g Drug %s',cellNumb,expDate,holdVol,drugCond);
uDirs = [uDirs; uDirs(1)];
uDirsPlot = repmat(uDirs,1,nReps);


%% Plot figure
% if ~showLess
%     if makeFig
        hF = figure;
%     else
%         hF = [];
%     end
%     
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

function [DSI] = getDSI(pSort,stimDir,prefDir)

ctMean = mean(pSort,2); %takes the mean of each row in the sorted ON or OFF responses = average response for each direction 
ctSortPlot = [pSort; pSort(1,:)]; %produces a matrix of the sorted responses ie is just pSortON or pSort OFF
ctMeanPlot = [ctMean; ctMean(1,:)]; %matrix of the means 

nReps = size(pSort,2); %defining number of repitions of directions based on number of columns in pSort 

uDirs = unique(stimDir);% list of directions %potentially problematic way of doing things;
% if stimDirs is shifted in a non-uniform way (e.g. resets to 0 after passing 360)
% then ctSort will now be incorrectly indexed

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
preSet_rlim = [0 1500];
prefSpikes = ctMean(prefIndx);
nullSpikes = ctMean(nullIndx);
titleStr = sprintf('Pref Dir %3.1f DSI %4.2f Vec Length %4.2f',prefDir,DSI,vecLength);
acqInfo = sprintf('Neuron %g Date %g Hold %g Drug %s',cellNumb,expDate,holdVol,drugCond);
uDirs = [uDirs; uDirs(1)];
uDirsPlot = repmat(uDirs,1,nReps);
end 