load cTableWThorz_FR.mat 
WTcTable = dTable; 
load cTableB2horz_FR.mat 
B2cTable = dTable; 

%DSI WT
data = [WTcTable.DSIon WTcTable.DSIoff];

figure('Name','DSI plot')
    tiledlayout(1,2);
    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('WT DSI')

%DSI B2
data = [B2cTable.DSIon B2cTable.DSIoff];

    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('B2 DSI')


%FR FOR WT
data = [WTcTable.prefFRon WTcTable.prefFRoff];

figure('Name','Max FR')
    tiledlayout(1,4);
    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Max Firing Rate (Hz)')
    ylim([0 100])
    title('WT Pref FR')

data = [WTcTable.nullFRon WTcTable.nullFRoff];
    
    nexttile 
    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Max Firing Rate (Hz)')
    ylim([0 100])
    title('WT Null FR')


%FR  B2
data = [B2cTable.prefFRon B2cTable.prefFRoff];
nexttile 
 plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Max Firing Rate (Hz)')
    ylim([0 100])
    title('B2 Pref FR')

    data = [B2cTable.nullFRon B2cTable.nullFRoff];
    
    nexttile 
    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Max Firing Rate (Hz)')
    ylim([0 100])
    title('B2 Null FR')

%     %Number of Spikes Wt
% data = [WTcTable.prefSpikesON WTcTable.prefSpikesOFF];
% 
% figure('Name','Number of Spikes')
%     tiledlayout(1,4);
%     nexttile 
% 
%     plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
%     'yLabel','# Spikes')
%     ylim([0 100])
%     title('WT Pref Spikes')
% 
%     data = [WTcTable.nullSpikesON WTcTable.nullSpikesOFF];
% 
%     nexttile 
%     plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
%     'yLabel','# Spikes')
%     ylim([0 100])
%     title('WT Null Spikes')
% 
%  %number spikes B2
% data = [B2cTable.prefSpikesON B2cTable.prefSpikesOFF];
% nexttile 
%  plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
%     'yLabel','# Spikes')
%     ylim([0 100])
%     title('B2 Pref Spikes')
% 
%     data = [B2cTable.nullSpikesON B2cTable.nullSpikesOFF];
% 
%     nexttile 
%     plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
%     'yLabel','# Spikes')
%     ylim([0 100])
%     title('B2 Null Spikes')
    


