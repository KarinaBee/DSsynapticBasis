%plot avg firing rate for all cells pref and null on and off and DSI FR

load cTableWTvert_FR.mat 
WTvTable = dTable; 

load cTableB2vert_FR.mat 
B2vTable = dTable;

load cTableB2horz_FR.mat 
B2hTable = dTable; 

load cTableWThorz_FR.mat
WThTable = dTable; 

%get values 
% WTvDSIon = WTvTable.DSIon;
% WTvDSIoff = WTvTable.DSIoff;
% WThDSIoff = WThTable.DSIoff;
% WThDSIon = WThTable.DSIon;
% B2hDSIon = B2hTable.DSIon;
% B2hDSIoff = B2hTable.DSIoff;
% B2vDSIoff = B2vTable.DSIoff;
% B2vDSIon = B2vTable.DSIon;
% WThPrefon = WThTable.prefFRon;
% WThNullon = WThTable.nullFRon;
% WThNulloff = WThTable.nullFRoff;
% WThPrefoff = WThTable.prefFRoff;
% 
% B2hPrefon = B2hTable.prefFRon;
% B2hNullon = B2hTable.nullFRon;
% B2hNulloff = B2hTable.nullFRoff;
% B2hPrefoff = B2hTable.prefFRoff;
% 
% WTvPrefon = WTvTable.prefFRon;
% WTvNullon = WTvTable.nullFRon;
% WTvNulloff = WTvTable.nullFRoff;
% WTvPrefoff = WTvTable.prefFRoff;
% 
% B2vPrefon = B2vTable.prefFRon;
% B2vNullon = B2vTable.nullFRon;
% B2vNulloff = B2vTable.nullFRoff;
% B2vPrefoff = B2vTable.prefFRoff;

%FR DSI WT Vert
data = [WTvTable.DSIon WTvTable.DSIoff];

figure('Name','DSI FR plot')
    tiledlayout(1,2);
    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('WT VERT FR DSI')

%FR DSI B2 Vert

data = [B2vTable.DSIon B2vTable.DSIoff];

    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('B2 VERT FR DSI')

%FR DSI WT Horz

data = [WThTable.DSIon WThTable.DSIoff];

    figure('Name','DSI FR plot')
    tiledlayout(1,2);
    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('WT HORZ FR DSI')
 
%FR DSI B2 Horz   

data = [B2hTable.DSIon B2hTable.DSIoff];
    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('B2 HORZ FR DSI')

%plotting pairwise FR for pref and null

%WT horz ON
    FRpref = WThTable.prefFRon;
    FRnull = WThTable.nullFRon;
   
    data = [FRpref FRnull];

    figure('Name','Max FR')
    tiledlayout(2,2);
    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(FRpref)
            plot([FRpref(k) FRnull(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
            xlim([0 3])
            ylim([0 100])
            ylabel('Max FR (Hz)')
            title('Max FR ON WT Horz')
            hold on
        end

 % WT HORZ OFF
    FRpref = WThTable.prefFRoff;
    FRnull = WThTable.nullFRoff;
   
data = [FRpref FRnull];
nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(FRpref)
            plot([FRpref(k) FRnull(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
            xlim([0 3])
            ylim([0 100])
            ylabel('Max FR (Hz)')
            title('Max FR OFF WT Horz')
            hold on
        end

% B2 HORZ ON
    FRpref = B2hTable.prefFRon;
    FRnull = B2hTable.nullFRon;
   
data = [FRpref FRnull];
 nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(FRpref)
            plot([FRpref(k) FRnull(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
            xlim([0 3])
            ylim([0 100])
            ylabel('Max FR (Hz)')
            title('Max FR ON B2 Horz')
            hold on
        end

 % B2 HORZ OFF
    FRpref = B2hTable.prefFRoff;
    FRnull = B2hTable.nullFRoff;
   
data = [FRpref FRnull];
 nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(FRpref)
            plot([FRpref(k) FRnull(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
            xlim([0 3])
            ylim([0 100])
            ylabel('Max FR (Hz)')
            title('Max FR OFF B2 Horz')
            hold on
        end

        
%WT VERT ON
    FRpref = WTvTable.prefFRon;
    FRnull = WTvTable.nullFRon;
   
    data = [FRpref FRnull];

    figure('Name','Max FR')
    tiledlayout(2,2);
    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(FRpref)
            plot([FRpref(k) FRnull(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
            xlim([0 3])
            ylim([0 200])
            ylabel('Max FR (Hz)')
            title('Max FR ON WT Vert')
            hold on
        end

 % WT HORZ OFF
    FRpref = WTvTable.prefFRoff;
    FRnull = WTvTable.nullFRoff;
   
    data = [FRpref FRnull];
 nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(FRpref)
            plot([FRpref(k) FRnull(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
            xlim([0 3])
            ylim([0 200])
            ylabel('Max FR (Hz)')
            title('Max FR OFF WT Vert')
            hold on
        end

% B2 HORZ ON
    FRpref = B2vTable.prefFRon;
    FRnull = B2vTable.nullFRoff;
   
    data = [FRpref FRnull];
 nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(FRpref)
            plot([FRpref(k) FRnull(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
            xlim([0 3])
            ylim([0 200])
            ylabel('Max FR (Hz)')
            title('Max FR ON B2 Vert')
            hold on
        end

        % B2 HORZ OFF
    FRpref = B2vTable.prefFRoff;
    FRnull = B2vTable.nullFRoff;
   
 data = [FRpref FRnull];
 nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(FRpref)
            plot([FRpref(k) FRnull(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
            xlim([0 3])
            ylim([0 200])
            ylabel('Max FR (Hz)')
            title('Max FR OFF B2 Vert')
            hold on
        end


% %plotting pairwise number of spikes for pref and null
% 
% %WT horz ON
%     FRpref = WThTable.prefSpikesON;
%     FRnull = WThTable.nullSpikesON;
% 
%     data = [FRpref FRnull];
% 
%     figure('Name','Max Spikes')
%     tiledlayout(2,2);
%     nexttile 
% 
%     plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
%         for k = 1 :size(FRpref)
%             plot([FRpref(k) FRnull(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
%             xlim([0 3])
%             ylim([0 200])
%             ylabel('Max FR (Hz)')
%             title('Max Spikes ON WT Horz')
%             hold on
%         end
% 
%  % WT HORZ OFF
%     FRpref = WThTable.prefSpikesOFF;
%     FRnull = WThTable.nullSpikesOFF;
% 
%     data = [FRpref FRnull];
%  nexttile 
% 
%     plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
%         for k = 1 :size(FRpref)
%             plot([FRpref(k) FRnull(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
%             xlim([0 3])
%             ylim([0 200])
%             ylabel('Max FR (Hz)')
%             title('Max Spikes OFF WT Horz')
%             hold on
%         end
% 
% % B2 HORZ ON
%     FRpref = B2hTable.prefSpikesON;
%     FRnull = B2hTable.nullSpikesON;
% 
%     data = [FRpref FRnull];
%  nexttile 
% 
%     plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
%         for k = 1 :size(FRpref)
%             plot([FRpref(k) FRnull(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
%             xlim([0 3])
%             ylim([0 200])
%             ylabel('Max FR (Hz)')
%             title('Max Spikes ON B2 Horz')
%             hold on
%         end
% 
%         % B2 HORZ OFF
%  FRpref = B2hTable.prefSpikesOFF;
%     FRnull = B2hTable.nullSpikesOFF;
% 
%     data = [FRpref FRnull];
%  nexttile 
% 
%     plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
%         for k = 1 :size(FRpref)
%             plot([FRpref(k) FRnull(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
%             xlim([0 3])
%             ylim([0 200])
%             ylabel('Max FR (Hz)')
%             title('Max Spikes OFF B2 Horz')
%             hold on
%         end
% 
%         %WT VERT ON
%     FRpref = WTvTable.prefSpikesON;
%     FRnull = WTvTable.nullSpikesON;
% 
%     data = [FRpref FRnull];
% 
%     figure('Name','Max Spikes')
%     tiledlayout(2,2);
%     nexttile 
% 
%     plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
%         for k = 1 :size(FRpref)
%             plot([FRpref(k) FRnull(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
%             xlim([0 3])
%             ylim([0 200])
%             ylabel('Max FR (Hz)')
%             title('Max Spikes ON WT Vert')
%             hold on
%         end
% 
%  % WT vert OFF
%     FRpref = WTvTable.prefSpikesOFF;
%     FRnull = WTvTable.nullSpikesOFF;
% 
%     data = [FRpref FRnull];
%  nexttile 
% 
%     plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
%         for k = 1 :size(FRpref)
%             plot([FRpref(k) FRnull(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
%             xlim([0 3])
%             ylim([0 200])
%             ylabel('Max FR (Hz)')
%             title('Max Spikes OFF WT Vert')
%             hold on
%         end
% 
% % B2 vert ON
%      FRpref = B2vTable.prefSpikesON;
%     FRnull = B2vTable.nullSpikesON;
% 
%     data = [FRpref FRnull];
%  nexttile 
% 
%     plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
%         for k = 1 :size(FRpref)
%             plot([FRpref(k) FRnull(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
%             xlim([0 3])
%             ylim([0 200])
%             ylabel('Max FR (Hz)')
%             title('Max Spikes ON B2 Vert')
%             hold on
%         end
% 
%         % B2 vert OFF
%          FRpref = B2vTable.prefSpikesOFF;
%     FRnull = B2vTable.nullSpikesOFF;
% 
%     data = [FRpref FRnull];
%  nexttile 
% 
%     plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
%         for k = 1 :size(FRpref)
%             plot([FRpref(k) FRnull(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');
%             xlim([0 3])
%             ylim([0 200])
%             ylabel('Max FR (Hz)')
%             title('Max spikes OFF B2 Vert')
%             hold on
%         end


% % Define the window bins
% windowBins = (0:0.1:1.4); % 15 bins of size 0.1 seconds
% 
% vWTavgFRonPD = mean(cell2mat(WTvTable.avgPDfiringRateON));
% vWTavgFRonND = mean(cell2mat(WTvTable.avgNDfiringRateON));
% 
% figure;
% hold on;
% bar(windowBins, vWTavgFRonPD, 'FaceAlpha', 0.5, 'DisplayName', 'Preferred Direction');
% bar(windowBins, vWTavgFRonND, 'FaceAlpha', 0.5, 'DisplayName', 'Null Direction');
% 
% % Customize the plot
% xlabel('Time (s)');
% ylabel('Firing Rate (Hz)');
% ylim([0 125])
% title('Firing Rates WT Vert ON');
% legend('show');
% 
% hold off;
% 
% vWTavgFRoffPD = mean(cell2mat(WTvTable.avgPDfiringRateOFF));
% vWTavgFRoffND = mean(cell2mat(WTvTable.avgNDfiringRateOFF));
% 
% figure;
% hold on;
% bar(windowBins, vWTavgFRoffPD, 'FaceAlpha', 0.5, 'DisplayName', 'Preferred Direction');
% bar(windowBins, vWTavgFRoffND, 'FaceAlpha', 0.5, 'DisplayName', 'Null Direction');
% 
% % Customize the plot
% xlabel('Time (s)');
% ylabel('Firing Rate (Hz)');
% ylim([0 125])
% title('Firing Rates WT Vert OFF');
% legend('show');
% 
% hold off;
% 
% vB2avgFRonPD = mean(cell2mat(B2vTable.avgPDfiringRateON));
% vB2avgFRonND = mean(cell2mat(B2vTable.avgNDfiringRateON));
% 
% figure;
% hold on;
% bar(windowBins, vB2avgFRonPD, 'FaceAlpha', 0.5, 'DisplayName', 'Preferred Direction');
% bar(windowBins, vB2avgFRonND, 'FaceAlpha', 0.5, 'DisplayName', 'Null Direction');
% 
% % Customize the plot
% xlabel('Time (s)');
% ylabel('Firing Rate (Hz)');
% ylim([0 125])
% title('Firing Rates B2 Vert ON');
% legend('show');
% 
% hold off;
% 
% vB2avgFRoffPD = mean(cell2mat(B2vTable.avgPDfiringRateOFF));
% vB2avgFRoffND = mean(cell2mat(B2vTable.avgNDfiringRateOFF));
% 
% figure;
% hold on;
% bar(windowBins, vB2avgFRoffPD, 'FaceAlpha', 0.5, 'DisplayName', 'Preferred Direction');
% bar(windowBins, vB2avgFRoffND, 'FaceAlpha', 0.5, 'DisplayName', 'Null Direction');
% 
% % Customize the plot
% xlabel('Time (s)');
% ylabel('Firing Rate (Hz)');
% ylim([0 125])
% title('Firing Rates B2 Vert OFF');
% legend('show');
% 
% hold off;
% 
% hB2avgFRonPD = mean(cell2mat(B2hTable.avgPDfiringRateON));
% hB2avgFRonND = mean(cell2mat(B2hTable.avgNDfiringRateON));
% 
% figure;
% hold on;
% bar(windowBins, hB2avgFRonPD, 'FaceAlpha', 0.5, 'DisplayName', 'Preferred Direction');
% bar(windowBins, hB2avgFRonND, 'FaceAlpha', 0.5, 'DisplayName', 'Null Direction');
% 
% % Customize the plot
% xlabel('Time (s)');
% ylabel('Firing Rate (Hz)');
% ylim([0 125])
% title('Firing Rates B2 Horz ON');
% legend('show');
% 
% hold off;
% 
% hB2avgFRoffPD = mean(cell2mat(B2hTable.avgPDfiringRateOFF));
% hB2avgFRoffND = mean(cell2mat(B2hTable.avgNDfiringRateOFF));
% 
% figure;
% hold on;
% bar(windowBins, hB2avgFRoffPD, 'FaceAlpha', 0.5, 'DisplayName', 'Preferred Direction');
% bar(windowBins, hB2avgFRoffND, 'FaceAlpha', 0.5, 'DisplayName', 'Null Direction');
% 
% % Customize the plot
% xlabel('Time (s)');
% ylabel('Firing Rate (Hz)');
% ylim([0 125])
% title('Firing Rates B2 Horz OFF');
% legend('show');
% 
% hold off;
% 
% hWTavgFRonPD = mean(cell2mat(WThTable.avgPDfiringRateON));
% hWTavgFRonND = mean(cell2mat(WThTable.avgNDfiringRateON));
% 
% figure;
% hold on;
% bar(windowBins, hWTavgFRonPD, 'FaceAlpha', 0.5, 'DisplayName', 'Preferred Direction');
% bar(windowBins, hWTavgFRonND, 'FaceAlpha', 0.5, 'DisplayName', 'Null Direction');
% 
% % Customize the plot
% xlabel('Time (s)');
% ylabel('Firing Rate (Hz)');
% ylim([0 125])
% title('Firing Rates WT Horz ON');
% legend('show');
% 
% hold off;
% 
% hWTavgFRoffPD = mean(cell2mat(WThTable.avgPDfiringRateOFF));
% hWTavgFRoffND = mean(cell2mat(WThTable.avgNDfiringRateOFF));
% 
% figure;
% hold on;
% bar(windowBins, hWTavgFRoffPD, 'FaceAlpha', 0.5, 'DisplayName', 'Preferred Direction');
% bar(windowBins, hWTavgFRoffND, 'FaceAlpha', 0.5, 'DisplayName', 'Null Direction');
% 
% % Customize the plot
% xlabel('Time (s)');
% ylabel('Firing Rate (Hz)');
% ylim([0 125])
% title('Firing Rates WT Horz OFF');
% legend('show');
% 
% hold off;